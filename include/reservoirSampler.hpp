#pragma once

#include "utils.hpp"
#include "collection.hpp"

#include "logging.hpp"

#include <unordered_set>
#include <string>

using namespace std::chrono;

template <
uint32_t t_factorization_block_size = 64*1024,
uint32_t t_estimator_block_size = 16,
uint32_t t_down_size = 256
>
struct reservoirSampler{

public: 	
    static std::string container_type()
    {
        return std::to_string(t_estimator_block_size);
    }
    
    static std::string container_file_name(collection& col, uint64_t size_in_bytes)
    {
        auto size_in_mb = size_in_bytes / (1024 * 1024);
        return col.path + "/index/" + container_type() + ".sdsl";
    }
public:
	 void sample(collection& col, std::vector<uint64_t>& rs, bool rebuild,size_t size_in_bytes) {

			auto start_total = hrclock::now();

			sdsl::read_only_mapper<8> text(col.file_map[KEY_TEXT]);
            auto n = text.size();

			// (1) create frequency estimates
			// try to load the estimates instead of recomputing
		   	// uint32_t down_size = 256;
			auto rs_name = container_file_name(col,size_in_bytes) + "-RSample-" +std::to_string(t_down_size);
			fixed_hasher<t_estimator_block_size> rk;

			uint64_t rs_size = (text.size()-t_estimator_block_size+1)/t_down_size;
			if (! utils::file_exists(rs_name) || rebuild) {		
				auto start = hrclock::now();
				LOG(INFO) << "\t" << "Building Reservoir sample with downsize: " << t_down_size;
				//make a reservoir sampler and keep it in RAM
				uint64_t count = 0;
				uint64_t skip = 0;
				std::random_device rd;
				std::mt19937 gen(rd());
				std::uniform_real_distribution<double> dis(0, 1);

				//double r = dis(gen);
				double w = std::exp(std::log(dis(gen))/rs_size);
				double s = std::floor(std::log(dis(gen))/std::log(1-w));
			
				for(size_t i=0;i<text.size();i++) {
					auto sym = text[i];
					auto hash = rk.update(sym);

					if(i < t_estimator_block_size-1) continue;
					else {
						if(count < rs_size) 
							rs.push_back(hash);
						else {
							skip++;
						// std::cout.precision(20);
						// std::cout << "\t" << "w=" << w << " s=" << s;	
						//	std::uniform_int_distribution<uint64_t> dis(0, count);
						//	uint64_t pos = dis(gen);
						//	if(pos < rs_size)
						//		rs[pos] = hash;
							if(s + 1 <= text.size()-t_estimator_block_size) {
								if(skip == s + 1) {
									// LOG(INFO) << "\t" << "here: " << count;
									rs[1+std::floor(rs_size*dis(gen))] = hash;
									w *= std::exp(std::log(dis(gen))/rs_size);
									s = std::floor(std::log(dis(gen))/std::log(1-w));
									skip = 0;
								}
							}
							else
								break; 	
						}
						count++;
					}						
				}
				auto stop = hrclock::now();
				LOG(INFO) << "\t" << "Reservoir sampling time = " << duration_cast<milliseconds>(stop-start).count() / 1000.0f << " sec";
				LOG(INFO) << "\t" << "Store reservoir sample to file " << rs_name;
				// sdsl::store_to_file(rs,rs_name);

				auto rs_file = sdsl::write_out_buffer<64>::create(rs_name);
	            {
	            	auto itr = rs.begin();
	            	while(itr != rs.end()) {
						std::copy(itr, itr+1, std::back_inserter(rs_file));
						itr++;
					}
	            }
			} else {
				LOG(INFO) << "\t" << "Load reservoir sample from file " << rs_name;
				// sdsl::load_from_file(rs,sketch_name);
				sdsl::read_only_mapper<64> rs_file(rs_name);
				auto itr = rs_file.begin();
				while(itr != rs_file.end()) {
					rs.push_back(*itr);
					itr++;
				}
			}
			
			LOG(INFO) << "\t" << "Reservoir sample blocks = " << rs.size(); 	
			LOG(INFO) << "\t" << "Reservoir sample size = " << rs.size()*8/(1024*1024) << " MiB";
	}
};
