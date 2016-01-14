#pragma once

#include "utils.hpp"
#include "collection.hpp"

#include "rlz_store_static.hpp"

template <class t_dictionary_creation_strategy,
          class t_dictionary_pruning_strategy,
          class t_dictionary_index,
          uint32_t t_factorization_block_size,
          bool t_search_local_block_context,
          class t_factor_selection_strategy,
          class t_factor_coder,
          class t_block_map>
class rlz_store_static<t_dictionary_creation_strategy,
                       t_dictionary_pruning_strategy,
                       t_dictionary_index,
                       t_factorization_block_size,
                       t_search_local_block_context,
                       t_factor_selection_strategy,
                       t_factor_coder,
                       t_block_map>::builder {
public:
    using dictionary_creation_strategy = t_dictionary_creation_strategy;
    using dictionary_pruning_strategy = t_dictionary_pruning_strategy;
    using dictionary_index_type = t_dictionary_index;
    using factor_selection_strategy = t_factor_selection_strategy;
    using factor_encoder = t_factor_coder;
    using factorization_strategy = factorizor<t_factorization_block_size,t_search_local_block_context, dictionary_index, factor_selection_strategy, factor_encoder>;
    using block_map = t_block_map;
    enum { block_size = t_factorization_block_size };
    enum { search_local_block_context = t_search_local_block_context };
public:
    builder& set_rebuild(bool r)
    {
        rebuild = r;
        return *this;
    };
    builder& set_threads(uint8_t nt)
    {
        num_threads = nt;
        return *this;
    };
    builder& set_dict_size(uint64_t ds)
    {
        dict_size_bytes = ds;
        return *this;
    };
    builder& set_pruned_dict_size(uint64_t ds)
    {
        pruned_dict_size_bytes = ds;
        return *this;
    };

    rlz_store_static build_or_load(collection& col) const
    {
        auto start = hrclock::now();

        // (1) create dictionary based on parametrized
        // dictionary creation strategy if necessary
        LOG(INFO) << "Create dictionary (" << dictionary_creation_strategy::type() << ")";
        dictionary_creation_strategy::create(col, rebuild, dict_size_bytes);
        LOG(INFO) << "Dictionary hash before pruning '" << col.param_map[PARAM_DICT_HASH] << "'";

        // (2) prune the dictionary if necessary
        LOG(INFO) << "Prune dictionary with " << dictionary_pruning_strategy::type();
        dictionary_pruning_strategy::template prune<dictionary_index_type, factorization_strategy>(col,
                                                                                                   rebuild, pruned_dict_size_bytes, num_threads);
        LOG(INFO) << "Dictionary after pruning '" << col.param_map[PARAM_DICT_HASH] << "'";

        // (3) create factorized text using the dict
        auto factor_file_name = factorization_strategy::factor_file_name(col);
        if (rebuild || !utils::file_exists(factor_file_name)) {
            factorization_strategy::template parallel_factorize<factor_storage>(col, rebuild, num_threads);
        } else {
            LOG(INFO) << "Factorized text exists.";
            col.file_map[KEY_FACTORIZED_TEXT] = factor_file_name;
            col.file_map[KEY_SPEP] = factorization_strategy::spep_file_name(col);
            col.file_map[KEY_BLOCKOFFSETS] = factorization_strategy::boffsets_file_name(col);
            col.file_map[KEY_BLOCKFACTORS] = factorization_strategy::bfactors_file_name(col);
            col.file_map[KEY_BLOCKSPEP] = factorization_strategy::bspep_file_name(col);
        }

        auto stop = hrclock::now();
        LOG(INFO) << "RLZ construction complete. time = " << duration_cast<seconds>(stop - start).count() << " sec";

        return rlz_store_static(col);
    }

    rlz_store_static load(collection& col) const
    {
        /* make sure components exists and register them */

        /* (1) check dict */
        auto dict_file_name = dictionary_creation_strategy::file_name(col, dict_size_bytes);
        if (!utils::file_exists(dict_file_name)) {
            throw std::runtime_error("LOAD FAILED: Cannot find dictionary.");
        } else {
            col.file_map[KEY_DICT] = dict_file_name;
            col.compute_dict_hash();
        }

        /* (2) check factorized text */
        auto dict_hash = col.param_map[PARAM_DICT_HASH];
        auto factor_file_name = factorization_strategy::factor_file_name(col);
        if (!utils::file_exists(factor_file_name)) {
            throw std::runtime_error("LOAD FAILED: Cannot find factorized text.");
        } else {
            col.file_map[KEY_FACTORIZED_TEXT] = factor_file_name;
            col.file_map[KEY_SPEP] = factorization_strategy::spep_file_name(col);
            col.file_map[KEY_BLOCKOFFSETS] = factorization_strategy::boffsets_file_name(col);
            col.file_map[KEY_BLOCKFACTORS] = factorization_strategy::bfactors_file_name(col);
            col.file_map[KEY_BLOCKSPEP] = factorization_strategy::bspep_file_name(col);
        }

        /* load */
        return rlz_store_static(col);
    }

    template<class t_idx>
    rlz_store_static remap_dict_apriori(t_idx& old,collection& col) const
    {
        auto start = hrclock::now();
        LOG(INFO) << "Remap dictionary offsets";
        /* (0) make sure we use the block size etc. */
        static_assert(t_idx::block_size == t_factorization_block_size,
                      "different factorization block size");
        static_assert(t_idx::factor_coder_type::literal_threshold == factor_coder_type::literal_threshold,
                      "different factor literal coding threshold");
        static_assert(t_idx::search_local_block_context == search_local_block_context,
                      "different local search strategy");
              
        LOG(INFO) << "Load dict index";
        dictionary_index_type idx(col,false);

        auto segment_size = dictionary_creation_strategy::segment_size();
        auto page_size = dictionary_creation_strategy::page_size();
        auto mode = dictionary_creation_strategy::remap_mode();
        uint8_t segment_size_log2 = sdsl::bits::hi(segment_size);
        uint8_t page_size_log2 = sdsl::bits::hi(page_size);

        /* (6) create page table */
        auto segment_mapping_file = col.path + "/index/SEGMENT-" + old.m_dict_hash + ".sdsl";
        sdsl::int_vector<> segment_table;
        sdsl::int_vector<> segment_table_rev;
        if(utils::file_exists(segment_mapping_file)) {
            LOG(INFO) << "Load segment mapping";
            std::ifstream ifs(segment_mapping_file);
            segment_table.load(ifs);
            segment_table_rev.load(ifs);
        } else {
            /* (3) iterate over factors to score segments */
            LOG(INFO) << "Score segments";
            LOG(INFO) << "Segment size = " << segment_size;
            
            auto fitr = old.factors_begin();
            auto fend = old.factors_end();
            uint64_t num_segments = dict_size_bytes / segment_size;
            std::vector<double> segment_score(num_segments);
            for(size_t i=0;i<segment_score.size();i++) segment_score[i] = 0;
            while (fitr != fend) {
                const auto& f = *fitr;
                if(!f.is_literal) {      
                    size_t m = f.ep-f.sp+1;
                    for(size_t i=0;i<m;i++) {
                        auto offset = idx.translate_offset(f.sp+i,f.len);
                        auto segment_id = offset / segment_size;
                        segment_score[segment_id] += 1.0f/(double)m;
                    }
                }
                ++fitr;
            }
            /* (4) sort scored segments */
            LOG(INFO) << "Sort segments by score";
            struct segment_info {
                double score;
                uint64_t id;
                segment_info(double s,uint64_t _id) : score(s), id(_id) {};
                bool operator<(segment_info& si) {
                    return si.score < score;
                } 
            };
            /* (5) sort scored segments */
            std::vector<segment_info> scored_segments;
            for(size_t i=0;i<segment_score.size();i++) {
                scored_segments.emplace_back(segment_score[i],i);
            }
            std::sort(scored_segments.begin(),scored_segments.end());
            
            LOG(INFO) << "Create page mapping";
            segment_table.resize(scored_segments.size());
            segment_table_rev.resize(scored_segments.size());
            for(size_t i=0;i<scored_segments.size();i++) {
                segment_table[i] = scored_segments[i].id;
                segment_table_rev[scored_segments[i].id] = i;
            }
            LOG(INFO) << "Store page mapping";
            std::ofstream ofs(segment_mapping_file);
            segment_table.serialize(ofs);
            segment_table_rev.serialize(ofs);
        }
        
        /* (7) create new dictionary */
        LOG(INFO) << "Create new dictionary";
        auto new_dict_file_name = dictionary_creation_strategy::file_name(col, dict_size_bytes);
        if(!utils::file_exists(new_dict_file_name)) {
            /* (1) get the old dict */
            const auto& existing_dict = old.dict;
            auto existing_dict_hash = old.m_dict_hash;
            col.file_map[KEY_DICT] = new_dict_file_name;
            auto new_dict = sdsl::write_out_buffer<8>::create(col.file_map[KEY_DICT]);
            /* (4) write the reordered dictionary */
            for(size_t i=0;i<segment_table.size();i++) {
                auto old_segment_id = segment_table[i];
                auto old_segment_itr = existing_dict.begin() + (old_segment_id*segment_size);
                for(size_t j=0;j<segment_size;j++) {
                    new_dict.push_back(*old_segment_itr++);
                }
            }
            new_dict.push_back(0); // EOD symbol
        } else {
            col.file_map[KEY_DICT] = new_dict_file_name;
        }
        col.compute_dict_hash();
        col.param_map[PARAM_DICT_HASH] = col.param_map[PARAM_DICT_HASH] + dictionary_creation_strategy::remap_mode_str();
        
        /* (8) map stuff */
        LOG(INFO) << "Remap factor offsets";
        {
            auto fitr = old.factors_begin();
            auto fend = old.factors_end();
            
            block_factor_data bfd(t_factorization_block_size);
            col.file_map[KEY_FACTORIZED_TEXT] = factorization_strategy::factor_file_name(col);
            col.file_map[KEY_SPEP] = factorization_strategy::spep_file_name(col);
            col.file_map[KEY_BLOCKOFFSETS] = factorization_strategy::boffsets_file_name(col);
            col.file_map[KEY_BLOCKFACTORS] = factorization_strategy::bfactors_file_name(col);
            col.file_map[KEY_BLOCKSPEP] = factorization_strategy::bspep_file_name(col);
            {
                auto factor_buf = sdsl::write_out_buffer<1>::create(col.file_map[KEY_FACTORIZED_TEXT]);
                auto spep_buf = sdsl::write_out_buffer<32>::create(col.file_map[KEY_SPEP]);
                auto block_offsets = sdsl::write_out_buffer<0>::create(col.file_map[KEY_BLOCKOFFSETS]);
                auto block_factors = sdsl::write_out_buffer<0>::create(col.file_map[KEY_BLOCKFACTORS]);
                auto block_spep = sdsl::write_out_buffer<0>::create(col.file_map[KEY_BLOCKSPEP]);
                bit_ostream<sdsl::int_vector_mapper<1> > factor_stream(factor_buf);
            
                t_factor_coder coder;
                auto num_blocks = old.block_map.num_blocks();
                auto num_blocks1p = (uint64_t)(num_blocks * 0.01);
    
                auto num_pages = dict_size_bytes / page_size;
                std::vector<uint8_t> preferred_pages(num_pages);
                std::vector<factor_data> cur_block_factors;
                while (fitr != fend) {
                    const auto& f = *fitr;
                    cur_block_factors.push_back(f);
                    if (fitr.last_factor_in_block()) {
                        /* process block previous block*/
                        process_block_factors(coder,
                            cur_block_factors,
                            bfd,
                            preferred_pages,
                            idx,
                            segment_table_rev,
                            mode,
                            page_size_log2,
                            segment_size_log2);
                
                        /* encode block  */
                        block_offsets.push_back(factor_stream.tellp());
                        block_factors.push_back(bfd.num_factors);
                        block_spep.push_back(spep_buf.size());
                        coder.encode_block(factor_stream,spep_buf,bfd);
                        bfd.reset();
                        cur_block_factors.clear();
                        if ((fitr.block_id + 1) % num_blocks1p == 0) {
                            LOG(INFO) << "\tEncoded " << 100 * (fitr.block_id + 1) / num_blocks << "% ("
                                    << (fitr.block_id + 1) << "/" << num_blocks << ")";
                        }
                    }
                    ++fitr;
                }
            }
        }
              
        auto stop = hrclock::now();
        LOG(INFO) << "RLZ factor remap complete. time = " << duration_cast<seconds>(stop - start).count() << " sec";
        return rlz_store_static(col);
    }
    
private:
    template<class t_coder,class t_idx>
    void process_block_factors(t_coder& coder,
        std::vector<factor_data>& cur_block_factors,
        block_factor_data& bfd,
        std::vector<uint8_t>& preferred_pages,
        t_idx& idx,
        sdsl::int_vector<>& segment_table,
        remap_mode_t mode,
        size_t page_size_log2,
        uint8_t segment_size_log2) const {
        for(size_t k=0;k<preferred_pages.size();k++) preferred_pages[k] = 0;
        if(mode == remap_mode_t::reorder_preferred_hapax) {
            /* add the hapaxes to preferred as we need them anyway */
            for(const auto& factor : cur_block_factors) {
                if(!factor.is_literal && factor.ep==factor.sp) {
                    auto cur_offset = idx.translate_offset(factor.sp,factor.len);
                    auto cur_segment = cur_offset >> segment_size_log2;
                    auto mapped_segment = segment_table[cur_segment];
                    auto mapped_page = mapped_segment >> page_size_log2;
                    preferred_pages[mapped_page] = 1;
                }
            }
        }
        
        uint32_t syms_in_block = 0;
        for(const auto& factor : cur_block_factors) {
            if(factor.is_literal) {
                bfd.add_factor(coder,factor.literal_ptr,0,factor.len,0,0);
            } else { 
                size_t m = factor.ep-factor.sp+1;
                bool segment_found = false;
                uint64_t in_page_offset = 0;
                uint64_t new_segment_id = 0;
                uint64_t min_segment_id = 99999999;
                uint64_t min_offset = 0;
                
                /* (1) map all offsets to segments */
                for(size_t i=0;i<m;i++) {
                    auto cur_offset = idx.translate_offset(factor.sp+i,factor.len);
                    auto cur_segment = cur_offset >> segment_size_log2;
                    auto mapped_segment = segment_table[cur_segment];
                    auto mapped_page = mapped_segment >> page_size_log2;
                    if(mode != remap_mode_t::only_reorder) {
                        if( preferred_pages[mapped_page] == 1 ) {
                            segment_found = true;
                            in_page_offset = cur_offset&sdsl::bits::lo_set[segment_size_log2];
                            new_segment_id = mapped_segment;
                            break;
                        }
                    }
                    if(mapped_segment < min_segment_id) {
                        min_segment_id = mapped_segment;
                        min_offset = cur_offset;
                    }
                }
            
                /* (2) if not found use smallest mapped segment */
                if(!segment_found) {
                    in_page_offset = min_offset&sdsl::bits::lo_set[segment_size_log2];
                    new_segment_id = min_segment_id;
                    auto new_page = new_segment_id >> page_size_log2;
                    if(mode != remap_mode_t::only_reorder) preferred_pages[new_page] = 1;
                }
                auto new_offset = (new_segment_id << segment_size_log2) + in_page_offset;
                bfd.add_factor(coder,factor.literal_ptr,new_offset,factor.len,factor.sp,factor.ep);
            }
            syms_in_block += factor.len;
        }  
    }
    
private:
    bool rebuild = false;
    uint32_t num_threads = 1;
    uint64_t dict_size_bytes = 0;
    uint64_t pruned_dict_size_bytes = 0;
};
