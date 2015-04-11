
#include "chunk_freq_estimator.hpp"
#include "sdsl/int_vector_mapper.hpp"

#include <cstdio>
#include <fstream>

#include <ctype.h>

template<class t_vec,class t_est>
void test_sketch(const t_vec& text,t_est& cfe)
{
	std::cout << "chunk size = " << cfe.chunk_size << std::endl;
	std::cout << "sketch size = " << cfe.sketch.sketch.size_in_bytes() << std::endl;
	std::cout << "d = " <<  cfe.sketch.sketch.d << std::endl;
	std::cout << "w = " <<  cfe.sketch.sketch.w << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    for(const auto& sym : text) {
        cfe.update(sym);
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto time_spent = stop-start;

    auto topk = cfe.topk();
    size_t j=1;
    for(const auto& item : topk) {
    	if(item.first == 1) break;
    	std::cout << std::setw(4) <<j<<" item = " << std::setw(10) << item.second.start << " ('";
    	auto itr = text.begin()+item.second.start;
    	auto end = itr + cfe.chunk_size;
    	while(itr != end) {
    		auto sym = *itr;
    		if(isspace(sym)) std::cout << " ";
    		else if(isprint(sym)) std::cout << sym;
    		else std::cout << "?";
    		++itr;
    	}
    	std::cout << "')";
        std::cout << " hash=" << std::setw(10) << item.second.hash << " efreq=" << std::setw(10) << item.first << std::endl;
        j++;
    }

    std::cout << "time in sec = " << std::chrono::duration_cast<std::chrono::milliseconds>(time_spent).count() / 1000.0f << std::endl;
}

int main(int argc,const char* argv[])
{
    if(argc!=2) return -1;
    auto file_name = argv[1];
    sdsl::read_only_mapper<8> text(file_name,true);

    {
		const uint32_t chunk_size = 8;
	    using sketch_type = count_min_sketch<std::ratio<1, 100000>,std::ratio<1, 10>>;
	    chunk_freq_estimator<chunk_size,1000,sketch_type> cfe;
	    test_sketch(text,cfe);
	}
    {
		const uint32_t chunk_size = 16;
	    using sketch_type = count_min_sketch<std::ratio<1, 100000>,std::ratio<1, 10>>;
	    chunk_freq_estimator<chunk_size,1000,sketch_type> cfe;
	    test_sketch(text,cfe);
	}
    {
		const uint32_t chunk_size = 32;
	    using sketch_type = count_min_sketch<std::ratio<1, 100000>,std::ratio<1, 10>>;
	    chunk_freq_estimator<chunk_size,1000,sketch_type> cfe;
	    test_sketch(text,cfe);
	}
    {
		const uint32_t chunk_size = 64;
	    using sketch_type = count_min_sketch<std::ratio<1, 100000>,std::ratio<1, 10>>;
	    chunk_freq_estimator<chunk_size,1000,sketch_type> cfe;
	    test_sketch(text,cfe);
	}
    {
		const uint32_t chunk_size = 128;
	    using sketch_type = count_min_sketch<std::ratio<1, 100000>,std::ratio<1, 10>>;
	    chunk_freq_estimator<chunk_size,1000,sketch_type> cfe;
	    test_sketch(text,cfe);
	}
    {
		const uint32_t chunk_size = 256;
	    using sketch_type = count_min_sketch<std::ratio<1, 100000>,std::ratio<1, 10>>;
	    chunk_freq_estimator<chunk_size,1000,sketch_type> cfe;
	    test_sketch(text,cfe);
	}
    {
		const uint32_t chunk_size = 512;
	    using sketch_type = count_min_sketch<std::ratio<1, 100000>,std::ratio<1, 10>>;
	    chunk_freq_estimator<chunk_size,1000,sketch_type> cfe;
	    test_sketch(text,cfe);
	}
    {
		const uint32_t chunk_size = 1024;
	    using sketch_type = count_min_sketch<std::ratio<1, 100000>,std::ratio<1, 10>>;
	    chunk_freq_estimator<chunk_size,1000,sketch_type> cfe;
	    test_sketch(text,cfe);
	}


    return 0;
}
