#pragma once

#include "count_min_sketch.hpp"

struct chunk_info {
	uint64_t start;
	uint64_t hash;
	bool operator==(const chunk_info& ci) const {
		return hash == ci.hash;
	}
	bool operator<(const chunk_info& ci) const {
		return start < ci.start;
	}
};

namespace std
{
    template<>
    struct hash<chunk_info>
    {
        typedef chunk_info argument_type;
        typedef std::size_t result_type;
         result_type operator()(argument_type const& s) const
        {
            return s.hash;
        }
    };
}

template<
uint32_t t_chunk_size = 256,
uint32_t t_k = 1000,
class t_sketch = count_min_sketch<>
>
struct chunk_freq_estimator {
public:
	using topk_sketch = count_min_sketch_topk<chunk_info,t_sketch>;
private:
	std::array<uint8_t,t_chunk_size> m_window; 
	uint64_t m_cur_offset = 0;
	uint64_t m_cur_KR_hash = 0;
	uint64_t num_chars = std::numeric_limits<uint8_t>::max() + 1ULL;
	uint64_t prime = 2147483647ULL; // largest prime <2^31
	uint64_t nk;
	topk_sketch m_freq_sketch = topk_sketch(t_k);
	std::queue<uint8_t> cur_chunk;
private:
	void compute_nk() {
		nk = 1;
		for(size_t i=0;i<t_chunk_size;i++) {
			nk = (nk*num_chars)%prime;
		}
	}
	chunk_info update_KR_hash(uint64_t tail,uint64_t front) {
		// /* (1) scale up */
		m_cur_KR_hash = (m_cur_KR_hash * num_chars)%prime;
		// /* (2) add last sym */
		m_cur_KR_hash = (m_cur_KR_hash + front)%prime;
		// /* (3) substract tail sym */
		m_cur_KR_hash = (m_cur_KR_hash + (prime - ((nk*tail)%prime)))%prime;
		return chunk_info{m_cur_offset-t_chunk_size,m_cur_KR_hash};
	}
public:
	uint32_t chunk_size = t_chunk_size;
	const topk_sketch& sketch = m_freq_sketch;
	chunk_freq_estimator() {
		compute_nk();
	}
	void update(uint8_t sym) {
		m_cur_offset++;
		cur_chunk.push(sym);
		if(m_cur_offset <= t_chunk_size) {
			m_cur_KR_hash = (m_cur_KR_hash*num_chars+sym)%prime;
			if(m_cur_offset == t_chunk_size) {
				chunk_info ci{0,m_cur_KR_hash};
				m_freq_sketch.update(ci);
			}
		} else {
			uint8_t tail = cur_chunk.front();
			cur_chunk.pop();
			auto ci = update_KR_hash(tail,sym);
			m_freq_sketch.update(ci);
		}
	}

	auto topk() -> decltype(m_freq_sketch.topk()) const  {
		return m_freq_sketch.topk();
	}
};