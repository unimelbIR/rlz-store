#pragma once

#include "utils.hpp"
#include "collection.hpp"

template <uint32_t t_block_size>
class init_dict_regularSampling {
public:
    static std::string type()
    {
        return "init_dict_regularSampling-"
               + std::to_string(t_block_size);
    }

    static std::string file_name(collection& col, uint64_t size_in_bytes)
    {
        auto size_in_mb = size_in_bytes / (1024 * 1024);
        return col.path + "/index/" + type() + "-" + std::to_string(size_in_mb) + ".sdsl";
    }

public:
    static void create(collection& col, bool rebuild, size_t size_in_bytes)
    {
        // const uint32_t block_size = t_block_size;
        uint64_t budget_bytes = size_in_bytes;
        uint64_t budget_mb = budget_bytes / (1024 * 1024);
        // check if we store it already and load it
        auto fname = file_name(col, size_in_bytes);
        col.file_map[KEY_DICT] = fname;
        if (!utils::file_exists(fname) || rebuild) { // construct
            auto start_total = hrclock::now();
            LOG(INFO) << "\tCreate dictionary with budget " << budget_mb << " MiB";
            // memory map the text and iterate over it
            sdsl::read_only_mapper<8> text(col.file_map[KEY_TEXT]);
            auto num_samples = budget_bytes / t_block_size;
            LOG(INFO) << "\tDictionary samples = " << num_samples;
            auto n = text.size();
            size_t sample_step = n / num_samples;   
            size_t sample_step_adjusted = sample_step / t_block_size * t_block_size;
            size_t num_samples_adjusted = n / sample_step_adjusted; //may contain more samples
            LOG(INFO) << "\tSample steps = " << sample_step;
            LOG(INFO) << "\tAdjusted sample steps = " << sample_step_adjusted;
            LOG(INFO) << "\tAdjusted dictionary samples  = " << num_samples_adjusted;

            std::vector<uint64_t> picked_blocks;
            std::unordered_set<uint64_t> segs;
            segs.max_load_factor(0.1);
            std::hash<std::string> hash_fn; 
            for(size_t i = 0; i < num_samples_adjusted; i++) {//steps
                uint64_t step_pos = i*sample_step_adjusted;
                for(size_t j = 0;j < sample_step_adjusted; j = j + t_block_size) {//blocks //assuming impossible to find all the same!!!
                    //impose hash to remove duplicates and kind increase diversitiy, a faster operation
                    auto start = text.begin()+step_pos + j;
                    std::string seg(start,start+t_block_size);              
                    size_t seg_hash = hash_fn(seg);
                    if(segs.find(seg_hash) != segs.end()) {
                        LOG(INFO) << "\t" << "same best block encountered!"; 
                        continue;
                    } else {
                        picked_blocks.push_back(step_pos + j);
                        segs.emplace(seg_hash);
                        break;
                    }
                }
                if(picked_blocks.size() >= num_samples) break; //breakout if dict is filled since adjusted is bigger
            }
           
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
            LOG(INFO) << "\t" << type() + " Total time = " << duration_cast<milliseconds>(end_total - start_total).count() / 1000.0f << " sec";
        } else {
            LOG(INFO) << "\t"
                      << "Dictionary exists at '" << fname << "'";
        }
        // compute a hash of the dict so we don't reconstruct things
        // later when we don't have to.
        col.compute_dict_hash();
    }
};
