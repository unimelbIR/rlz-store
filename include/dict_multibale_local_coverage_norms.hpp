#pragma once

#include "utils.hpp"
#include "collection.hpp"

#include "logging.hpp"

#include "count_min_sketch.hpp"
#include "chunk_freq_estimator.hpp"

#include <unordered_set>

using namespace std::chrono;

template <
uint32_t t_block_size = 1024,
uint32_t t_estimator_block_size = 16,
uint32_t t_down_size = 512,
class t_norm = std::ratio<1,2>
>
class dict_multibale_local_coverage_norms{
public:
    static std::string type()
    {
        return "dict_multibale_local_coverage_norms-" +std::to_string(t_block_size)+"-"+ std::to_string(t_estimator_block_size);
    }
    static uint32_t adjusted_down_size(collection& col, uint64_t size_in_bytes)
    {
		sdsl::read_only_mapper<8> text(col.file_map[KEY_TEXT]);
		auto thres = (text.size()/size_in_bytes)/2;
        // return (thres >= t_down_size? thres : t_down_size);
		return 512;   
 	}
 	
    static std::string container_type()
    {
        return std::to_string(t_estimator_block_size);
    }
    static std::string dict_file_name(collection& col, uint64_t size_in_bytes, int c_type, int isCombined = 0)
    {
        auto size_in_mb = size_in_bytes / (1024 * 1024);
		std::string ctype = "-rw" + std::to_string(c_type);
		return col.path + "/index/"  + col.text + "-" + type() + "-" + std::to_string(size_in_mb) + "-" + std::to_string(t_norm::num) + "-" + std::to_string(t_norm::den) + + "-" + std::to_string(adjusted_down_size(col, size_in_bytes))+".sdsl"+ctype+"-"+std::to_string(isCombined);
    }
    static std::string container_file_name(collection& col, uint64_t size_in_bytes)
    {
        auto size_in_mb = size_in_bytes / (1024 * 1024);
        return col.path + "/index/"  + col.text + "-" + container_type() + ".sdsl";
    }
public:
	static void create(collection& col, bool rebuild,size_t size_in_bytes, int c_type, std::unordered_set<uint64_t> &history_mers, int isCombined = 0) {
		uint32_t budget_bytes = size_in_bytes;
		uint32_t budget_mb = size_in_bytes / (1024 * 1024);
		// uint32_t num_blocks_required = budget_bytes / t_block_size;

        // check if we store it already and load it
        auto fname = dict_file_name(col, size_in_bytes, c_type, isCombined);
        col.file_map[KEY_DICT] = fname;
		auto down_size = adjusted_down_size(col, size_in_bytes);
		if (! utils::file_exists(fname) || rebuild ) {  // construct
			auto start_total = hrclock::now();
			LOG(INFO) << "\t" << "Create dictionary with budget " << budget_mb << " MiB";
			LOG(INFO) << "\t" << "Block size = " << t_block_size; 
			// LOG(INFO) << "\t" << "Num blocks = " << num_blocks_required; 

			sdsl::read_only_mapper<8> text(col.file_map[KEY_TEXT]);       
            auto n = text.size();
            size_t num_samples = budget_bytes / t_block_size; //hopefully much smaller than the adjusted
            size_t scale =  n / budget_bytes; //hopefully much smaller than the adjusted, may not be divisible, can fix later
            size_t sample_step = scale * t_block_size;
		    int thres = 1; 
		    if(scale >= 8*1024) thres = 16;
		    else if(scale >= 1024 && scale < 8*1024) thres = 8;
		    else if(scale >= 512 && scale < 1024) thres = 4; //
		    else thres = 2; //minimum saving half

            size_t sample_step_adjusted = sample_step / thres / t_block_size * t_block_size; //make a tempate para later
            size_t num_samples_adjusted = n / sample_step_adjusted; //may contain more samples

            // size_t sample_step = n / num_samples;   
            // size_t sample_step_adjusted = sample_step / t_block_size * t_block_size;
            // size_t num_samples_adjusted = n / sample_step_adjusted; //may contain more samples

            //fix sampling every 1mb pick 1kb until dictionary is filled

            LOG(INFO) << "\tDictionary samples to be picked = " << num_samples;
            LOG(INFO) << "\tText epoch size = " << sample_step;
            LOG(INFO) << "\tAdjusted epoch size = " << sample_step_adjusted;
            LOG(INFO) << "\tAdjusted dictionary samples in the text = " << num_samples_adjusted;

			// (1) create frequency estimates
			// try to load the estimates instead of recomputing
		   	// uint32_t down_size = 256;
			auto rs_name = container_file_name(col,size_in_bytes) + "-RSample-" +std::to_string(down_size);
			fixed_hasher<t_estimator_block_size> rk;

			uint64_t rs_size = text.size()/down_size;
			std::vector<uint64_t> rs; //filter out frequency less than 64
			if (! utils::file_exists(rs_name) || rebuild) {		
				auto start = hrclock::now();
				LOG(INFO) << "\t" << "Building Reservoir sample with downsize: " << down_size;
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
			// LOG(INFO) << "\t" << "Reservoir sample block hashes = " << rs; 
			LOG(INFO) << "\t" << "Reservoir sample blocks = " << rs.size(); 	
			LOG(INFO) << "\t" << "Reservoir sample size = " << rs.size()*8/(1024*1024) << " MiB";
			//build exact counts of sampled elements
			LOG(INFO) << "\t" << "Calculating exact frequencies of small rolling blocks...";
			std::unordered_map<uint64_t,uint32_t> mers_counts;
			mers_counts.max_load_factor(0.1);
			for(uint64_t s : rs) {
				mers_counts[s]++;	
			}
			rs.clear(); //might be able to do it in place!!!!

			//auto itr = rs.begin();
			// while(rs.begin() != rs.end()) {
			// 	mers_counts[*rs.begin()]++;	
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
			LOG(INFO) << "\t" << "Useful kept small blocks no. = " << mers_counts.size(); 	

			//first pass getting densest steps!
			// std::vector<std::pair <uint32_t,uint64_t>> steps;
			LOG(INFO) << "\t" << "First pass: getting random steps..."; 
			auto start = hrclock::now();

			std::vector<uint32_t> step_indices;
			//try ordered max cov from front to back, proven to be good for small datasets
			for (size_t i = 0; i < num_samples_adjusted; i++) {
				step_indices.push_back(i);
			}
			// //try randomly ordered max cov
			std::random_shuffle(step_indices.begin(), step_indices.end());
			
		    auto stop = hrclock::now();
		    LOG(INFO) << "\t" << "1st pass runtime = " << duration_cast<milliseconds>(stop-start).count() / 1000.0f << " sec";

			// 2nd pass: process max coverage using the sorted order by density
			std::unordered_set<uint64_t> step_mers; //can be prefilled?
			step_mers.max_load_factor(0.1);

			//prefill step_mers from the GOV2 bales history_mers?.sdsl hard code this file name, to be replaced by Matt later
			// std::string bale_file = "history_mers_BALE";
			// std::string bale_no = "5";//number to be changed, windows should be auto as well.

			// std::string histroy_file = col.path + "/index/" + bale_file + bale_no; //to be manually changed

			// LOG(INFO) << "\t" << "Load history_mers from file " << histroy_file; 
			// sdsl::load_from_file(rs,sketch_name);
			// sdsl::read_only_mapper<64> history_mapper(histroy_file);
			if(!history_mers.empty()) {
				auto itr = history_mers.begin();
				while(itr != history_mers.end()) {
					step_mers.emplace(*itr);
					itr++;
				}
				LOG(INFO) << "\t" << "Prefill history done with total mers size " << step_mers.size(); 
			} else {
				LOG(INFO) << "\t" << "No mers to prefill"; 
			}
			
			//process max cov
			std::vector<uint64_t> picked_blocks;
			LOG(INFO) << "\t" << "Second pass: perform ordered max coverage..."; 
			start = hrclock::now();
			double norm = (double)t_norm::num/t_norm::den;
			LOG(INFO) << "\t" << "Computing norm = " << norm; 
			for(const auto& i : step_indices) {//steps
				double sum_weights_max = std::numeric_limits<double>::min();	
				uint64_t step_pos = i*sample_step_adjusted;
				uint64_t best_block_no = step_pos;
				std::unordered_set<uint64_t> best_local_mers;

				for(size_t j=0;j<sample_step_adjusted;j = j+t_block_size) {//blocks 
					std::unordered_set<uint64_t> local_mers;
					local_mers.max_load_factor(0.1);
					double sum_weights_current = 0;

					//computational expensive place
					for(size_t k=0;k<t_block_size;k++) {//bytes?k=k+2?
						auto sym = text[step_pos+j+k];
						auto hash = rk.update(sym);

						if(k < t_estimator_block_size-1) continue;

						if(local_mers.find(hash) == local_mers.end() && step_mers.find(hash) == step_mers.end()) //continues rolling
						{//expensive checking		
							if(mers_counts.find(hash) != mers_counts.end()) {
								local_mers.emplace(hash);
								auto freq = mers_counts[hash];
								//compute norms
								sum_weights_current += std::pow(freq,norm); //L0.5
							}			
						}
					}
					//sum_weights_current = std::sqrt(sum_weights_current);
					if(norm > 0)
						sum_weights_current = std::pow(sum_weights_current,1/norm);
					if(sum_weights_current >= sum_weights_max) //consider tight breaking may get marginal gain, maynot be worth it.
					{
						sum_weights_max = sum_weights_current;
						best_block_no = step_pos + j; 
						best_local_mers = std::move(local_mers);
					}
				}
			
				picked_blocks.push_back(best_block_no);
				// LOG(INFO) << "\t" << "Blocks picked: " << picked_blocks.size(); 
				if(picked_blocks.size() >= num_samples) break; //breakout if dict is filled since adjusted is bigger
				step_mers.insert(best_local_mers.begin(), best_local_mers.end());
				// //erase selected blocks
				// for(auto hash : best_local_mers) {
				// 	mers_counts.erase(hash);
				// }
		    }   
			LOG(INFO) << "\t" << "Blocks size to check = " << step_mers.size(); 
		    step_mers.clear(); //save mem
		    step_indices.clear();//
		    // useful_blocks.clear();
		    mers_counts.clear();	
		    //sort picked blocks
			std::sort(picked_blocks.begin(), picked_blocks.end());
			stop = hrclock::now();
			LOG(INFO) << "\t" << "Picked blocks = " << picked_blocks; 
			LOG(INFO) << "\t" << "2nd pass runtime = " << duration_cast<milliseconds>(stop-start).count() / 1000.0f << " sec";

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
