#pragma once

#include "utils.hpp"
#include "collection.hpp"

#include "logging.hpp"

#include <unordered_set>
#include <string>

using namespace std::chrono;

template <
uint32_t t_estimator_block_size = 1024
>
class deduplicator{
public:
 	static std::string type()
    {
        return "dict_deduplicator-"+ std::to_string(t_estimator_block_size);
    }

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
	static void create(collection& col, bool rebuild,size_t size_in_bytes) {

			auto start_total = hrclock::now();
			sdsl::read_only_mapper<8> text(col.file_map[KEY_TEXT]);
            auto n = text.size();
			
			auto d_name = container_file_name(col,size_in_bytes) + "-DistinctDisjoint-";
			uint64_t d_size = text.size()/t_estimator_block_size;

			std::unordered_set<uint64_t> distinctBlocks;
			distinctBlocks.max_load_factor(0.1);

			//change back to roll hashes for later comparison
			if (! utils::file_exists(d_name) || rebuild) {		
				auto start = hrclock::now();
				LOG(INFO) << "\t" << "First pass: counting disjoint duplicates...";
				uint64_t count = 0;
				for(size_t i=0;i<text.size()-t_estimator_block_size+1;i=i+t_estimator_block_size) {
                    fixed_hasher<t_estimator_block_size> rk;	
                    uint64_t hash = -1;

                    for(size_t j=0;j<t_estimator_block_size;j++) {
                    	auto sym = text[i+j];
						hash = rk.update(sym);
					}

                    if(distinctBlocks.find(hash) == distinctBlocks.end()) {
                    	distinctBlocks.emplace(hash);
                    } 
                    else {
                    	count++;
                    	LOG(INFO) << "\t" << "Disjoint Duplicated blocks: " << count;
                    }
	                LOG(INFO) << "\t" << "Completed: " << (double)i*100/text.size() << "\%"; 
				}
				auto stop = hrclock::now();
				LOG(INFO) << "\t" << "Disjoint Deduplication time = " << duration_cast<milliseconds>(stop-start).count() / 1000.0f << " sec";
            	LOG(INFO) << "\t" << "Disjoint Duplication in kb: " << count*t_estimator_block_size/1024;
            	LOG(INFO) << "\t" << "Disjoint Duplication percentage: " << (double)count*100/d_size << "\%";

				LOG(INFO) << "\t" << "Store distinct disjoint block hashes to file " << d_name;

				// sdsl::store_to_file(rs,rs_name);

				auto rs_file = sdsl::write_out_buffer<64>::create(d_name);
	            {
	            	auto itr = distinctBlocks.begin();																																																										
	            	while(itr != distinctBlocks.end()) {
						std::copy(itr, itr+1, std::back_inserter(rs_file));
						itr++;
					}
	            }
			} else {
				LOG(INFO) << "\t" << "Load distinct disjoint blocks from file " << d_name;
				// sdsl::load_from_file(rs,sketch_name);
				sdsl::read_only_mapper<64> d_file(d_name);
				auto itr = d_file.begin();
				while(itr != d_file.end()) {
					distinctBlocks.emplace(*itr);
					itr++;
				}
			}
			
			//2nd pass
			/*LOG(INFO) << "\t" << "Second pass: counting rolling duplicates by skipping...";
			uint64_t count = 0;
			fixed_hasher<t_estimator_block_size> rk;
			for(size_t i=0;i<text.size();i++) {
				auto hash = rk.update(sym);

				if(i < t_estimator_block_size-1) continue;

                if(distinctBlocks.find(hash) == distinctBlocks.end()) {//look for the next rolling block
                	continue;
                } 
                else {//skip to the next non-overlap item
                	count++;
                	LOG(INFO) << "\t" << "Disjoint Duplicated blocks: " << count;
                }
                LOG(INFO) << "\t" << "Completed: " << (double)i*100/text.size() << "\%"; 
			}
			auto stop = hrclock::now();
			LOG(INFO) << "\t" << "Disjoint Deduplication time = " << duration_cast<milliseconds>(stop-start).count() / 1000.0f << " sec";
        	LOG(INFO) << "\t" << "Disjoint Duplication in kb: " << count*t_estimator_block_size/1024;
        	LOG(INFO) << "\t" << "Disjoint Duplication percentage: " << (double)count*100/d_size << "\%";*/
	}
};
