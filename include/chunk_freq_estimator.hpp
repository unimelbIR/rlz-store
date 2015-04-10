#pragma once

struct chunk_info {
	uint64_t start;
	uint64_t hash;
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
private:
	std::array<uint8_t,t_chunk_size> m_window; 
	uint64_t m_cur_offset = 0;
	uint64_t m_cur_KR_hash = 4711;
	uint64_t num_chars = std::numeric_limits<uint8_t>::max() + 1ULL;
	uint64_t prime = 2147483647ULL; // largest prime <2^32
	uint64_t nk;
	count_min_sketch_topk<chunk_info,t_sketch> m_freq_sketch(t_k);
	std::queue<uint8_t> cur_chunk;
private:
	void compute_nk() {
		nk = 1;
		for(size_t i=0;i<t_chunk_size;i++) {
			nk = (nk*num_chars)&prime;
		}
	}
	chunk_info update_KR_hash(uint8_t tail,uint8_t front) {
		/* (1) scale up */
		m_cur_KR_hash = (m_cur_KR_hash * num_chars)&prime;
		/* (2) add last sym */
		m_cur_KR_hash = (m_cur_KR_hash + front)&prime;
		/* (3) substract tail sym */
		auto sub = prime - ((nk*tail)&prime);
		m_cur_KR_hash = (m_cur_KR_hash + sub)&prime;
		return chunk_info{m_cur_offset,m_cur_KR_hash};
	}
public:
	template<class t_itr>
	chunk_freq_estimator(t_itr itr) {
		/* init the hash code */
		auto end = itr + t_chunk_size;
		while(itr != end) {
			auto sym = *itr;
			cur_chunk.push(sym);
			m_cur_KR_hash = (m_cur_KR_hash*num_chars+sym)&prime;
			itr++;
			m_cur_offset++;
		}
	}
	void update(uint8_t sym) {
		cur_chunk.push(sym);
		m_cur_offset++;
		uint8_t tail = cur_chunk.front();
		cur_chunk.pop();
		auto ci = update_KR_hash(tail,sym);
		m_freq_sketch.update(ci);
	}

	std::vector<chunk_info> topk() {
		return m_freq_sketch.topk();
	}
}