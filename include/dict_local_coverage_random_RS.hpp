#pragma once

#include "utils.hpp"
#include "collection.hpp"

#include "logging.hpp"

#include "count_min_sketch.hpp"
#include "chunk_freq_estimator.hpp"

#include <unordered_set>

using namespace std::chrono;
enum ACCESS_TYPE : int {
    SEQ,
    RAND
};

template <uint32_t t_block_size = 1024,
          uint32_t t_estimator_block_size = 16,
          uint32_t t_down_size = 256,
          ACCESS_TYPE t_method = SEQ>
class dict_local_coverage_random_RS {
public:
    static std::string container_type()
    {
        return std::to_string(t_estimator_block_size);
    }
    static std::string type()
    {
        return "dict_local_coverage_random_RS-" + std::to_string(t_method) + "-" + std::to_string(t_block_size) + "-" + std::to_string(t_estimator_block_size);
    }
    static std::string dict_file_name(collection& col, uint64_t size_in_bytes)
    {
        auto size_in_mb = size_in_bytes / (1024 * 1024);
        return col.path + "/index/" + type() + "-" + std::to_string(size_in_mb) + ".sdsl";
    }
    static std::string container_file_name(collection& col, uint64_t size_in_bytes)
    {
        auto size_in_mb = size_in_bytes / (1024 * 1024);
        return col.path + "/index/" + container_type() + ".sdsl";
    }

public:
    static void create(collection& col, bool rebuild, size_t size_in_bytes)
    {
        uint32_t budget_bytes = size_in_bytes;
        uint32_t budget_mb = size_in_bytes / (1024 * 1024);
        uint32_t num_blocks_required = (budget_bytes / t_estimator_block_size) + 1;

        // check if we store it already and load it
        auto fname = dict_file_name(col, size_in_bytes);
        col.file_map[KEY_DICT] = fname;
        if (!utils::file_exists(fname) || rebuild) { // construct
            auto start_total = hrclock::now();
            LOG(INFO) << "\t"
                      << "Create dictionary with budget " << budget_mb << " MiB";
            LOG(INFO) << "\t"
                      << "Block size = " << t_block_size;
            LOG(INFO) << "\t"
                      << "Num blocks = " << num_blocks_required;

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

            // (1) Filter by threshold using reservoir sampling
            // try to load the estimates instead of recomputing
            uint32_t down_size = t_down_size;
            auto rs_name = container_file_name(col, size_in_bytes) + "-RSample-" + std::to_string(down_size);
            fixed_hasher<t_estimator_block_size> rk;

            uint64_t rs_size = (text.size() - t_estimator_block_size) / down_size;
            std::vector<uint64_t> rs; //filter out frequency less than 64
            if (!utils::file_exists(rs_name)) { //to add || rebuild
                auto start = hrclock::now();
                LOG(INFO) << "\t"
                          << "Building Reservoir sample";
                //make a reservoir sampler and keep it in RAM
                uint64_t count = 0;
                uint64_t skip = 0;
                std::random_device rd;
                std::mt19937 gen(rd());
                std::uniform_real_distribution<double> dis(0, 1);

                //double r = dis(gen);
                double w = std::exp(std::log(dis(gen)) / rs_size);
                double s = std::floor(std::log(dis(gen)) / std::log(1 - w)); //may use floor which is more stable?

                for (size_t i = 0; i < text.size(); i++) {
                    auto sym = text[i];
                    auto hash = rk.update(sym);

                    if (i < t_estimator_block_size - 1)
                        continue;
                    else {
                        if (count < rs_size)
                            rs.push_back(hash);
                        else {
                            // std::uniform_int_distribution<uint64_t> dis(0, count);
                            // uint64_t pos = dis(gen);
                            // if(pos < rs_size)
                            // 	rs[pos] = hash;
                            skip++;
                            if (s + 1 <= text.size() - t_estimator_block_size) {
                                if (skip == s + 1) {
                                    // LOG(INFO) << "\t" << "w = " << w << "s = " << s;
                                    // LOG(INFO) << "\t" << "here: " << count;
                                    rs[1 + std::floor(rs_size * dis(gen))] = hash;
                                    w *= std::exp(std::log(dis(gen)) / rs_size);
                                    s = std::floor(std::log(dis(gen)) / std::log(1 - w));
                                    skip = 0;
                                }
                            } else
                                break;
                        }
                        count++;
                    }
                }
                auto stop = hrclock::now();
                LOG(INFO) << "\t"
                          << "Reservoir sampling time = " << duration_cast<milliseconds>(stop - start).count() / 1000.0f << " sec";
                LOG(INFO) << "\t"
                          << "Store reservoir sample to file " << rs_name;
                // sdsl::store_to_file(rs,rs_name);

                auto rs_file = sdsl::write_out_buffer<64>::create(rs_name);
                {
                    auto itr = rs.begin();
                    while (itr != rs.end()) {
                        std::copy(itr, itr + 1, std::back_inserter(rs_file));
                        itr++;
                    }
                }
            } else {
                LOG(INFO) << "\t"
                          << "Load reservoir sample from file " << rs_name;
                // sdsl::load_from_file(rs,sketch_name);
                sdsl::read_only_mapper<64> rs_file(rs_name);
                auto itr = rs_file.begin();
                while (itr != rs_file.end()) {
                    rs.push_back(*itr);
                    itr++;
                }
            }

            LOG(INFO) << "\t"
                      << "Reservoir sample blocks = " << rs.size();
            LOG(INFO) << "\t"
                      << "Reservoir sample size = " << rs.size() * 8 / (1024 * 1024) << " MiB";

            // (2) make a useful blocks hash, clear RA as you go to save space
            std::sort(rs.begin(), rs.end());
            auto last = std::unique(rs.begin(), rs.end());
            rs.erase(last, rs.end());
            // LOG(INFO) << "\t" << "RS size = " << rs.size();
            std::unordered_set<uint64_t> useful_blocks;
            useful_blocks.max_load_factor(0.2);
            // auto itr = rs.begin();
            // while(itr != rs.end()) {
            // 	useful_blocks.emplace(*itr);
            // 	itr = rs.erase(itr);
            // 	// itr++;
            // }
            std::move(rs.begin(), rs.end(), std::inserter(useful_blocks, useful_blocks.end()));
            LOG(INFO) << "\t"
                      << "Useful kept small blocks no. = " << useful_blocks.size();

            //1st pass getting densest steps!
            // std::vector<std::pair <uint32_t,uint64_t>> steps;
            LOG(INFO) << "\t"
                      << "First pass: getting random steps...";
            auto start = hrclock::now();

            std::vector<uint32_t> step_indices;
            //try ordered max cov from front to back, proven to be good for small datasets
            for (size_t i = 0; i < num_samples_adjusted; i++) {
                step_indices.push_back(i);
            }
            // //try randomly ordered max cov
            if (t_method == RAND)
                std::random_shuffle(step_indices.begin(), step_indices.end());

            auto stop = hrclock::now();
            LOG(INFO) << "\t"
                      << "1st pass runtime = " << duration_cast<milliseconds>(stop - start).count() / 1000.0f << " sec";

            // 2nd pass: process max coverage using the sorted order by density
            std::unordered_set<uint64_t> step_blocks;
            step_blocks.max_load_factor(0.2);
            std::vector<uint64_t> picked_blocks;
            LOG(INFO) << "\t"
                      << "Second pass: perform ordered max coverage...";
            start = hrclock::now();

            for (const auto& i : step_indices) { //steps
                uint32_t sum_weights_max = std::numeric_limits<uint32_t>::min();
                uint32_t intersection_max = 0;
                uint64_t step_pos = i * sample_step_adjusted;
                uint64_t best_block_no = step_pos;
                std::unordered_set<uint64_t> best_local_blocks;

                for (size_t j = 0; j < sample_step_adjusted; j = j + t_block_size) { //blocks
                    std::unordered_set<uint64_t> local_blocks;
                    uint32_t old_pos = 0;
                    uint32_t intersection_count = 0;
                    uint32_t count = 0;
                    local_blocks.max_load_factor(0.1);
                    uint32_t sum_weights_current = 0;

                    for (size_t k = 0; k < t_block_size; k++) { //bytes
                        auto sym = text[step_pos + j + k];
                        auto hash = rk.update(sym);

                        if (k < t_estimator_block_size - 1)
                            continue;

                        if (local_blocks.find(hash) == local_blocks.end() && step_blocks.find(hash) == step_blocks.end() && useful_blocks.find(hash) != useful_blocks.end()) //continues rolling
                        { //expensive checking
                            count++;
                            if (count > 1 && k - t_estimator_block_size <= old_pos) {
                                intersection_count++;
                            }
                            old_pos = k;
                            sum_weights_current += 1;
                            local_blocks.emplace(hash);
                        }
                    }
                    // LOG(INFO) << "\t" << "intersection_count = " << intersection_count;
                    // the way: partially resolve tights!!! prefer more intersections
                    if (sum_weights_current > sum_weights_max || (sum_weights_current == sum_weights_max && intersection_count > intersection_max)) {
                        //						LOG(INFO) << "\t" << "Rolling tight!";
                        intersection_max = intersection_count;
                        sum_weights_max = sum_weights_current;
                        best_block_no = step_pos + j;
                        best_local_blocks = std::move(local_blocks);
                    }
                }
                //				 LOG(INFO) << "\t" << "Block " << best_block_no << "\t" << "max weight = " << sum_weights_max;
                //				 LOG(INFO) << "\t" << "Block " << best_block_no << "\t" << "max intersect = " << intersection_max << "\n";

                picked_blocks.push_back(best_block_no);
                //				LOG(INFO) << "\t" << "\n";
                // LOG(INFO) << "\t" << "Blocks picked: " << picked_blocks.size();
                if (picked_blocks.size() >= num_samples)
                    break; //breakout if dict is filled since adjusted is bigger
                step_blocks.insert(best_local_blocks.begin(), best_local_blocks.end());
            }
            LOG(INFO) << "\t"
                      << "Blocks size to check = " << step_blocks.size();
            step_blocks.clear(); //save mem
            step_indices.clear(); //
            useful_blocks.clear();

            //sort picked blocks
            std::sort(picked_blocks.begin(), picked_blocks.end());
            stop = hrclock::now();
            LOG(INFO) << "\t"
                      << "Picked blocks = " << picked_blocks;
            LOG(INFO) << "\t"
                      << "2nd pass runtime = " << duration_cast<milliseconds>(stop - start).count() / 1000.0f << " sec";

            // last pass: writing to dict
            LOG(INFO) << "\t"
                      << "Last: writing dictionary...";
            auto dict = sdsl::write_out_buffer<8>::create(col.file_map[KEY_DICT]);
            {
                sdsl::read_only_mapper<8> text(col.file_map[KEY_TEXT]);
                for (const auto& pb : picked_blocks) {
                    auto beg = text.begin() + pb;
                    auto end = beg + t_block_size;
                    std::copy(beg, end, std::back_inserter(dict));
                }
            }

            LOG(INFO) << "\t"
                      << "Final dictionary size = " << dict.size() / (1024 * 1024) << " MiB";
            dict.push_back(0); // zero terminate for SA construction
            auto end_total = hrclock::now();
            LOG(INFO) << "\n"
                      << "\t" << type() + " Total time = " << duration_cast<milliseconds>(end_total - start_total).count() / 1000.0f << " sec";
        } else {
            LOG(INFO) << "\t"
                      << "Dictionary exists at '" << fname << "'";
        }
        // compute a hash of the dict so we don't reconstruct things
        // later when we don't have to.
        col.compute_dict_hash();
    }
};
