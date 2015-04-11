#pragma once

#include <ratio>
#include <cmath>
#include <random>
#include <queue>
#include <unordered_set>

#include <sdsl/int_vector.hpp>

struct hash_params {
	uint64_t a;
	uint64_t b;
};

template<
class t_epsilon = std::ratio<1, 20000>, // default settings ~ 2.5MB 
class t_delta = std::ratio<1, 1000>
>
struct count_min_sketch {
private:
	sdsl::int_vector<32> m_table;
	std::vector<hash_params> m_hash_params;
private:
	inline void pick_hash_params() {
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<uint64_t> param_dist(1, prime);
		for(size_t i=0;i<d;i++) {
			auto a = param_dist(gen); 
			auto b = param_dist(gen);
			m_hash_params.push_back({a,b});
		}
	}
	inline uint32_t compute_hash(uint64_t x,size_t row) const {
		uint64_t hash = (m_hash_params[row].a*x) + m_hash_params[row].b;
		hash = ((hash >> 31) + hash) & prime;
		return hash & w;
	}
public:
	double epsilon = (double) t_epsilon::num / (double) t_epsilon::den;
	double delta = (double) t_delta::num / (double) t_delta::den;
	uint64_t w_real = std::ceil(2.0 / epsilon);
	uint64_t w = (1ULL << (sdsl::bits::hi(w_real)+1ULL))-1ULL; // smallest power of 2 larger than w
	uint64_t d = std::ceil(std::log(1.0/delta)/std::log(2.0));
	uint64_t prime = 2147483647; // largest prime <2^32
	uint64_t total_count = 0;

	count_min_sketch() {
		m_table = sdsl::int_vector<32>(w*d);
		pick_hash_params();
	} 
	uint64_t update(uint64_t item,size_t count = 1) {
		total_count += count;
		uint64_t new_est = std::numeric_limits<uint64_t>::max();
		for(size_t i=0;i<d;i++) {
			auto row_offset = compute_hash(item,i);
			auto col_offset = w*i;
			uint64_t new_count = m_table[row_offset+col_offset] + count;
			m_table[row_offset+col_offset] = new_count;
			new_est = std::min(new_est,new_count);
		}
		return new_est;
	}
	uint64_t estimate(uint64_t item) const {
		uint64_t est = std::numeric_limits<uint64_t>::max();
		for(size_t i=0;i<d;i++) {
			auto row_offset = compute_hash(item,i);
			auto col_offset = w*i;
			uint64_t val = m_table[row_offset+col_offset];
			est = std::min(est,val);
		}
		return est;
	}
	uint64_t size_in_bytes() const {
		uint64_t bytes = 0;
		bytes += sdsl::size_in_bytes(m_table);
		bytes += d*sizeof(hash_params);
		return bytes;
	}
	double estimation_error() const {
		return epsilon * total_count;
	}
	double estimation_probability() const {
		return 1.0 - delta;
	}
};

template<
class T,
class t_cms = count_min_sketch<>
>
struct count_min_sketch_topk {
private:
	struct topk_item {
		T item;
		uint64_t estimate;
		topk_item(T& i,uint64_t e) : item(i), estimate(e) {};
		bool operator<(const topk_item& i) const {
			return estimate < i.estimate;
		}
		bool operator>(const topk_item& i) const {
			return estimate > i.estimate;
		}
		bool operator==(const T& i) {
			return item == i;
		}
	};
	t_cms m_sketch;
	std::vector<topk_item> m_topk_pq;
	std::unordered_set<T> m_topk_set;
	std::hash<T> hash_fn;
	uint64_t m_k;
private:
	inline void update_topk(T& item,uint64_t estimate) {
		if( m_topk_set.size() < m_k ) {
			auto itr = m_topk_set.find(item);
			if(itr == m_topk_set.end()) {
				m_topk_pq.emplace_back(item,estimate);
				m_topk_set.emplace(item);
				if(m_topk_pq.size() == m_k) {
					// create the heap after we have seen k items
					std::make_heap(m_topk_pq.begin(),m_topk_pq.end(),
								   std::greater<topk_item>());
				}
			}
		} else {
			auto itr = m_topk_set.find(item);
			if(itr != m_topk_set.end()) {
				// found the item. update estimate and ensure heap correctness
				auto pq_itr = m_topk_pq.rbegin();
				auto pq_end = m_topk_pq.rend();
				while(pq_itr != pq_end) {
					if(*pq_itr == item) {
						auto dist_from_end = std::distance(m_topk_pq.rbegin(),pq_itr);
						auto dist_from_front = m_topk_pq.size() - dist_from_end - 1;

						pq_itr->estimate = estimate;
						std::make_heap(m_topk_pq.begin()+dist_from_front,
							m_topk_pq.end(),std::greater<topk_item>());
						break;
					}
					++pq_itr;
				}
			} else {
				// check against smallest in pq
				if( estimate != 1 && m_topk_pq[0].estimate < estimate ) {
					// throw out the smallest and put the new one in
					std::pop_heap(m_topk_pq.begin(),m_topk_pq.end());
					m_topk_set.erase(m_topk_pq.back().item);
					m_topk_pq.pop_back();
					m_topk_pq.emplace_back(item,estimate);
					std::push_heap(m_topk_pq.begin(), m_topk_pq.end());
					m_topk_set.insert(item);
				}
			}
		}
	}
public:
	const t_cms& sketch = m_sketch;
	count_min_sketch_topk(size_t k) :m_k(k) {

	}
	void update(T& item,size_t count = 1) {
		uint64_t hash = hash_fn(item);
		auto estimate = m_sketch.update(hash,count);
		update_topk(item,estimate);
	}
	uint64_t estimate(T& item) const {
		uint64_t hash = hash_fn(item);
		return m_sketch.estimate(hash);
	}
	std::vector<std::pair<uint64_t,T>> topk() const {
		std::vector<std::pair<uint64_t,T>> cur_topk;
		for(const auto& pqi : m_topk_pq) {
			cur_topk.emplace_back(pqi.estimate,pqi.item);
		}
		std::sort(cur_topk.begin(),cur_topk.end(),std::greater<std::pair<uint64_t,T>>());
		return cur_topk;
	}
};