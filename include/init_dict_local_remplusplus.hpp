#pragma once

#include "utils.hpp"
#include "collection.hpp"

#include "logging.hpp"

#include "count_min_sketch.hpp"
#include "chunk_freq_estimator.hpp"
#include "reservoirSampler.hpp"

#include <unordered_set>

using namespace std::chrono;

template <
uint32_t t_block_size = 1024,
uint32_t t_estimator_block_size = 16,
uint32_t t_down_size = 256,
class t_norm = std::ratio<1,2>
>
class init_dict_local_remplusplus{
public:
    static std::string type()
    {
        return "init_dict_local_remplusplus-"+std::to_string(t_block_size)+"-"+ std::to_string(t_estimator_block_size);
    }
 	
    static std::string dict_file_name(collection& col, uint64_t size_in_bytes)
    {
        auto size_in_mb = size_in_bytes / (1024 * 1024);
		return col.path + "/index/" + type() + "-" + std::to_string(size_in_mb) + "-" + std::to_string(t_norm::num) + "-" + std::to_string(t_norm::den) + + "-" + std::to_string(adjusted_down_size(col, size_in_bytes))+".sdsl";
    }
    static uint32_t adjusted_down_size(collection& col, uint64_t size_in_bytes)
    {
		sdsl::read_only_mapper<8> text(col.file_map[KEY_TEXT]);
		//auto ratio = (text.size()/size_in_bytes)/2;
        //return (ratio >= t_down_size? t_down_size : ratio);
		return 256; //for 10gb
    }
   
public:
	static void create(collection& col, bool rebuild,size_t size_in_bytes) {
		uint32_t budget_bytes = size_in_bytes;
		uint32_t budget_mb = size_in_bytes / (1024 * 1024);
		uint32_t num_blocks_required = (budget_bytes / t_estimator_block_size) + 1;
		uint32_t down_size = adjusted_down_size(col, size_in_bytes);
        // check if we store it already and load it
        auto fname = dict_file_name(col, size_in_bytes);
        col.file_map[KEY_DICT] = fname;
		// auto down_size = adjusted_down_size(col, size_in_bytes);
		if (! utils::file_exists(fname) || rebuild ) {  // construct
			auto start_total = hrclock::now();
			LOG(INFO) << "\t" << "Create dictionary with budget " << budget_mb << " MiB";
			LOG(INFO) << "\t" << "Block size = " << t_block_size; 
			LOG(INFO) << "\t" << "Num blocks = " << num_blocks_required; 

			sdsl::read_only_mapper<8> text(col.file_map[KEY_TEXT]);
			auto num_samples = budget_bytes / t_block_size;
            LOG(INFO) << "\tDictionary samples = " << num_samples;
            auto n = text.size();
            size_t sample_step = n / num_samples;   
            size_t sample_step_adjusted = sample_step / t_block_size * t_block_size;
            size_t num_samples_adjusted = n / sample_step_adjusted; //may contain more samples
            fixed_hasher<t_estimator_block_size> rk;

            LOG(INFO) << "\tSample steps = " << sample_step;
            LOG(INFO) << "\tAdjusted sample steps = " << sample_step_adjusted;
            LOG(INFO) << "\tAdjusted dictionary samples  = " << num_samples_adjusted;

            //reservoir sampling
	        std::vector<uint64_t> rs;
            reservoirSampler<16,256> sampler;
	        sampler.sample(col, rs, rebuild, size_in_bytes);  

			//build exact counts of sampled elements
			LOG(INFO) << "\t" << "Calculating exact frequencies of small rolling blocks...";
			std::unordered_map<uint64_t,uint32_t> block_counts;
			block_counts.max_load_factor(0.1);
			// LOG(INFO) <<  (block_counts[12233] == 0);

			for(uint64_t s : rs) {
				block_counts[s]++;
			}
			rs.clear();
			//auto itr = rs.begin();
			// while(rs.begin() != rs.end()) {
			// 	block_counts[*rs.begin()]++;	
			// 	rs.erase(rs.begin());
			// }
			// LOG(INFO) << "\t" << "Generating distinct small rolling samples...";
			// //make a useful blocks hash, clear RA as you go to save space
			// std::sort(rs.begin(), rs.end());
			// auto last = std::unique(rs.begin(), rs.end());
	  //   		rs.erase(last, rs.end());
	  //   		// LOG(INFO) << "\t" << "RS size = " << rs.size(); 	
			// std::unordered_set<uint64_t> useful_blocks;		
			// useful_blocks.max_load_factor(0.2);
			
			// std::move(rs.begin(), rs.end(), std::inserter(useful_blocks, useful_blocks.end()));
			LOG(INFO) << "\t" << "Useful kept small blocks no. = " << block_counts.size();

			//REM++
          	std::hash<std::string> hash_fn;	
          	std::unordered_set<uint64_t> segs;
          	std::vector<uint64_t> picked_blocks;
          	auto start = hrclock::now();
			double norm = (double)t_norm::num/t_norm::den;
			LOG(INFO) << "\t" << "Computing norm = " << norm; 

			for (size_t i = 0; i < num_samples_adjusted; i++) {
				double sum_weights_max = std::numeric_limits<double>::min();	
				uint64_t step_pos = i*sample_step_adjusted;
				uint64_t best_block_no = step_pos;
				std::unordered_set<uint64_t> best_local_blocks;

				for(size_t j=0;j<sample_step_adjusted;j = j+t_block_size) {//blocks 
					std::unordered_set<uint64_t> local_blocks;
					double sum_weights_current = 0;
					auto start = text.begin()+step_pos + j;
                    std::string seg(start,start+t_block_size);				
                    size_t seg_hash = hash_fn(seg);
                    if(segs.find(seg_hash) != segs.end()) {
                    	LOG(INFO) << "\t" << "same best block encountered!"; 
                    	continue;
                    }
					for(size_t k=0;k<t_block_size;k++) {//bytes
						auto sym = text[step_pos+j+k];
						auto hash = rk.update(sym);

						if(k < t_estimator_block_size-1) continue;

						if(local_blocks.find(hash) == local_blocks.end() && block_counts.find(hash) != block_counts.end()) //continues rolling
						// if(local_blocks.find(hash) == local_blocks.end() && block_counts[hash] != 0) //continues rolling
						{//expensive checking					
							local_blocks.emplace(hash);
							auto freq = block_counts[hash];
							//compute norms
							sum_weights_current += std::pow(freq,norm); //L0.5
						}
					}
					//sum_weights_current = std::sqrt(sum_weights_current);
					if(norm > 0)
						sum_weights_current = std::pow(sum_weights_current,1/norm);
					if(sum_weights_current >= sum_weights_max)			
					{
						sum_weights_max = sum_weights_current;
						best_block_no = step_pos + j; 
					}
				}
				auto start = text.begin()+best_block_no;
                std::string seg(start,start+t_block_size);				
                size_t seg_hash = hash_fn(seg);
				segs.emplace(seg_hash);

				picked_blocks.push_back(best_block_no);
				// LOG(INFO) << "\t" << "Blocks picked: " << picked_blocks.size(); 
				if(picked_blocks.size() >= num_samples) break; //breakout if dict is filled since adjusted is bigger
		    }   

		    block_counts.clear();	
		    //sort picked blocks
			std::sort(picked_blocks.begin(), picked_blocks.end());
			auto stop = hrclock::now();
			LOG(INFO) << "\t" << "Picked blocks = " << picked_blocks; 
			LOG(INFO) << "\t" << "REM plusplus runtime = " << duration_cast<milliseconds>(stop-start).count() / 1000.0f << " sec";

		    // last pass: writing to dict
			LOG(INFO) << "\t" << "Last: writing dictionary..."; 
			auto dict = sdsl::write_out_buffer<8>::create(col.file_map[KEY_DICT]);
            {
			    sdsl::read_only_mapper<8> text(col.file_map[KEY_TEXT]);
			    for(const auto& pb : picked_blocks) {
				    auto beg = text.begin() + pb;
				    auto end = beg + t_block_size;
				    std::copy(beg,end,std::back_inserter(dict));
			    }
            }
	
			LOG(INFO) << "\t" << "Final dictionary size = " << dict.size()/(1024*1024) << " MiB"; 
			dict.push_back(0); // zero terminate for SA construction
			auto end_total = hrclock::now();
			LOG(INFO) << "\t" << type() << " Total time = " << duration_cast<milliseconds>(end_total-start_total).count() / 1000.0f << " sec";
		} else {
			LOG(INFO) << "\t" << "Dictionary exists at '" << fname << "'";
		}
		// compute a hash of the dict so we don't reconstruct things
		// later when we don't have to.
		col.compute_dict_hash();
	}
};
