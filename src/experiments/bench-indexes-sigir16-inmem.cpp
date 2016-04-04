#define ELPP_THREAD_SAFE
#define ELPP_STL_LOGGING

#include "utils.hpp"
#include "collection.hpp"
#include "rlz_utils.hpp"

#include "indexes.hpp"

#include "logging.hpp"
INITIALIZE_EASYLOGGINGPP

template<class t_idx,class t_dict_idx>
sdsl::int_vector<> create_or_load_segment_table(collection& col,t_idx& old,t_dict_idx& idx,size_t segment_size,size_t dict_size_bytes) {

    /* (6) create page table */
    auto segment_mapping_file = col.path + "/index/SEGMENT-" + old.m_dict_hash + ".sdsl";
    sdsl::int_vector<> segment_table;
    sdsl::int_vector<> segment_table_rev;
    if (utils::file_exists(segment_mapping_file)) {
        LOG(INFO) << "Load segment mapping";
        std::ifstream ifs(segment_mapping_file);
        segment_table.load(ifs);
        segment_table_rev.load(ifs);
    }
    else {
        /* (3) iterate over factors to score segments */
        LOG(INFO) << "Score segments";
        LOG(INFO) << "Segment size = " << segment_size;

        auto fitr = old.factors_begin();
        auto fend = old.factors_end();
        uint64_t num_segments = dict_size_bytes / segment_size;
        std::vector<double> segment_score(num_segments);
        for (size_t i = 0; i < segment_score.size(); i++)
            segment_score[i] = 0;
        while (fitr != fend) {
            const auto& f = *fitr;
            if (!f.is_literal) {
                size_t m = f.ep - f.sp + 1;
                for (size_t i = 0; i < m; i++) {
                    auto offset = idx.translate_offset(f.sp + i, f.len);
                    auto segment_id = offset / segment_size;
                    segment_score[segment_id] += 1.0f / (double)m;
                }
            }
            ++fitr;
        }
        /* (4) sort scored segments */
        LOG(INFO) << "Sort segments by score";
        struct segment_info {
            double score;
            uint64_t id;
            segment_info(double s, uint64_t _id)
                : score(s)
                , id(_id){};
            bool operator<(segment_info& si)
            {
                return si.score < score;
            }
        };
        /* (5) sort scored segments */
        std::vector<segment_info> scored_segments;
        for (size_t i = 0; i < segment_score.size(); i++) {
            scored_segments.emplace_back(segment_score[i], i);
        }
        std::sort(scored_segments.begin(), scored_segments.end());

        LOG(INFO) << "Create page mapping";
        segment_table.resize(scored_segments.size());
        segment_table_rev.resize(scored_segments.size());
        for (size_t i = 0; i < scored_segments.size(); i++) {
            segment_table[i] = scored_segments[i].id;
            segment_table_rev[scored_segments[i].id] = i;
        }
        LOG(INFO) << "Store page mapping";
        std::ofstream ofs(segment_mapping_file);
        segment_table.serialize(ofs);
        segment_table_rev.serialize(ofs);
    }
    return segment_table_rev;
}

template<class t_coder,class t_stream>
void compress_rlz_block_inmem(
    const char* name,
    t_coder& coder,
    size_t block_id,
    t_stream& os,
    std::vector<factor_data>& cur_block_factors,
    block_factor_data& bfd)
{
    /* reset data */
    os.seek(0);
    bfd.reset();
    
    /* parse data */
    for (const auto& factor : cur_block_factors) {
        if (factor.is_literal) {
            bfd.add_factor(coder, factor.literal_ptr, 0, factor.len, 0, 0, 0, 0);
        }
        else {
            bfd.add_factor(coder, factor.literal_ptr, factor.offset, factor.len, factor.sp, factor.ep, factor.sp1, factor.ep1);
        }
    }
    
    /* encode data  */
    coder.encode_block_output(os,bfd,block_id,name);
}

template<class t_coder_pv,class t_coder_zzz,class t_stream,class t_idx>
void compress_rlzp2a_block_inmem(
    t_coder_pv& coder_pv,
    t_coder_zzz& coder_zzz,
    size_t block_id,
    t_stream& os,
    std::vector<factor_data>& cur_block_factors,
    sdsl::int_vector<>& segment_table,
    uint64_t segment_size_log2,
    t_idx& idx,
    block_factor_data& bfd,
    block_factor_data& bfd_zzz)
{
    /* reset data */
    os.seek(0);
    bfd.reset();
    bfd_zzz.reset();
    
    /* parse data */
    for (const auto& factor : cur_block_factors) {
        if (factor.is_literal) {
            bfd.add_factor(coder_pv, factor.literal_ptr, 0, factor.len, 0, 0, 0, 0);
            bfd_zzz.add_factor(coder_zzz, factor.literal_ptr, 0, factor.len, 0, 0, 0, 0);
        }
        else {
            size_t m = factor.ep - factor.sp + 1;
            bool segment_found = false;
            uint64_t in_page_offset = 0;
            uint64_t new_segment_id = 0;
            uint64_t min_segment_id = 99999999;
            uint64_t min_offset = 0;

            /* (1) map all offsets to segments */
            for (size_t i = 0; i < m; i++) {
                auto cur_offset = idx.translate_offset(factor.sp + i, factor.len);
                auto cur_segment = cur_offset >> segment_size_log2;
                auto mapped_segment = segment_table[cur_segment];
                if (mapped_segment < min_segment_id) {
                    min_segment_id = mapped_segment;
                    min_offset = cur_offset;
                }
            }

            /* (2) if not found use smallest mapped segment */
            if (!segment_found) {
                in_page_offset = min_offset & sdsl::bits::lo_set[segment_size_log2];
                new_segment_id = min_segment_id;
            }
            auto new_offset = (new_segment_id << segment_size_log2) + in_page_offset;
            bfd.add_factor(coder_pv, factor.literal_ptr, new_offset, factor.len, factor.sp, factor.ep, factor.sp1, factor.ep1);
            bfd_zzz.add_factor(coder_zzz, factor.literal_ptr, new_offset, factor.len, factor.sp, factor.ep, factor.sp1, factor.ep1);
        }
    }
    
    /* encode data  */
    coder_pv.encode_block_output(os,bfd,block_id,"RLZP2A-PV");
    coder_zzz.encode_block_output(os,bfd_zzz,block_id,"RLZP2A-ZZZ");
}

template<class t_coder_pv,class t_coder_zzz,class t_stream,class t_idx>
void compress_rlzp2b_block_inmem(
    t_coder_pv& coder_pv,
    t_coder_zzz& coder_zzz,
    size_t block_id,
    t_stream& os,
    std::vector<factor_data>& cur_block_factors,
    sdsl::int_vector<>& segment_table,
    uint64_t segment_size_log2,
    uint64_t page_size_log2,
    std::vector<uint8_t>& preferred_pages,
    t_idx& idx,
    block_factor_data& bfd,
    block_factor_data& bfd_zzz)
{
    /* reset data */
    os.seek(0);
    bfd.reset();
    bfd_zzz.reset();
    for (size_t k = 0; k < preferred_pages.size(); k++) preferred_pages[k] = 0;
    
    /* parse data */
    for (const auto& factor : cur_block_factors) {
        if (factor.is_literal) {
            bfd.add_factor(coder_pv, factor.literal_ptr, 0, factor.len, 0, 0, 0, 0);
            bfd_zzz.add_factor(coder_zzz, factor.literal_ptr, 0, factor.len, 0, 0, 0, 0);
        }
        else {
            size_t m = factor.ep - factor.sp + 1;
            bool segment_found = false;
            uint64_t in_page_offset = 0;
            uint64_t new_segment_id = 0;
            uint64_t min_segment_id = 99999999;
            uint64_t min_offset = 0;

            /* (1) map all offsets to segments */
            for (size_t i = 0; i < m; i++) {
                auto cur_offset = idx.translate_offset(factor.sp + i, factor.len);
                auto cur_segment = cur_offset >> segment_size_log2;
                auto mapped_segment = segment_table[cur_segment];
                auto mapped_page =  (mapped_segment << segment_size_log2) >> page_size_log2;
                if (preferred_pages[mapped_page] == 1) {
                    segment_found = true;
                    in_page_offset = cur_offset & sdsl::bits::lo_set[segment_size_log2];
                    new_segment_id = mapped_segment;
                    break;
                }
                if (mapped_segment < min_segment_id) {
                    min_segment_id = mapped_segment;
                    min_offset = cur_offset;
                }
            }

            /* (2) if not found use smallest mapped segment */
            if (!segment_found) {
                in_page_offset = min_offset & sdsl::bits::lo_set[segment_size_log2];
                new_segment_id = min_segment_id;
                auto new_page = (new_segment_id << segment_size_log2) >> page_size_log2;
                preferred_pages[new_page] = 1;
            }
            auto new_offset = (new_segment_id << segment_size_log2) + in_page_offset;
            bfd.add_factor(coder_pv, factor.literal_ptr, new_offset, factor.len, factor.sp, factor.ep, factor.sp1, factor.ep1);
            bfd_zzz.add_factor(coder_zzz, factor.literal_ptr, new_offset, factor.len, factor.sp, factor.ep, factor.sp1, factor.ep1);
        }
    }
    
    /* encode data  */
    coder_pv.encode_block_output(os,bfd,block_id,"RLZP2B-PV");
    coder_zzz.encode_block_output(os,bfd_zzz,block_id,"RLZP2B-ZZZ");
}


template<class t_coder_pv,class t_coder_zzz,class t_stream,class t_idx>
void compress_rlzp2c_block_inmem(
    t_coder_pv& coder_pv,
    t_coder_zzz& coder_zzz,
    size_t block_id,
    t_stream& os,
    std::vector<factor_data>& cur_block_factors,
    sdsl::int_vector<>& segment_table,
    uint64_t segment_size_log2,
    uint64_t page_size_log2,
    std::vector<uint8_t>& preferred_pages,
    t_idx& idx,
    block_factor_data& bfd,
    block_factor_data& bfd_zzz)
{
    /* reset data */
    os.seek(0);
    bfd.reset();
    bfd_zzz.reset();
    for (size_t k = 0; k < preferred_pages.size(); k++) preferred_pages[k] = 0;
    
    /* add the hapaxes to preferred as we need them anyway */
    for (const auto& factor : cur_block_factors) {
        if (!factor.is_literal && factor.ep == factor.sp) {
            auto cur_offset = idx.translate_offset(factor.sp, factor.len);
            auto cur_segment = cur_offset >> segment_size_log2;
            auto mapped_segment = segment_table[cur_segment];
            auto mapped_page = (mapped_segment << segment_size_log2) >> page_size_log2;
            preferred_pages[mapped_page] = 1;
        }
    }
    
    /* parse data */
    for (const auto& factor : cur_block_factors) {
        if (factor.is_literal) {
            bfd.add_factor(coder_pv, factor.literal_ptr, 0, factor.len, 0, 0, 0,0);
            bfd_zzz.add_factor(coder_zzz, factor.literal_ptr, 0, factor.len, 0, 0, 0,0);
        }
        else {
            size_t m = factor.ep - factor.sp + 1;
            bool segment_found = false;
            uint64_t in_page_offset = 0;
            uint64_t new_segment_id = 0;
            uint64_t min_segment_id = 99999999;
            uint64_t min_offset = 0;

            /* (1) map all offsets to segments */
            for (size_t i = 0; i < m; i++) {
                auto cur_offset = idx.translate_offset(factor.sp + i, factor.len);
                auto cur_segment = cur_offset >> segment_size_log2;
                auto mapped_segment = segment_table[cur_segment];
                auto mapped_page =  (mapped_segment << segment_size_log2) >> page_size_log2;
                if (preferred_pages[mapped_page] == 1) {
                    segment_found = true;
                    in_page_offset = cur_offset & sdsl::bits::lo_set[segment_size_log2];
                    new_segment_id = mapped_segment;
                    break;
                }
                if (mapped_segment < min_segment_id) {
                    min_segment_id = mapped_segment;
                    min_offset = cur_offset;
                }
            }

            /* (2) if not found use smallest mapped segment */
            if (!segment_found) {
                in_page_offset = min_offset & sdsl::bits::lo_set[segment_size_log2];
                new_segment_id = min_segment_id;
                auto new_page = (new_segment_id << segment_size_log2) >> page_size_log2;
                preferred_pages[new_page] = 1;
            }
            auto new_offset = (new_segment_id << segment_size_log2) + in_page_offset;
            bfd.add_factor(coder_pv, factor.literal_ptr, new_offset, factor.len, factor.sp, factor.ep, factor.sp1, factor.ep1);
            bfd_zzz.add_factor(coder_zzz, factor.literal_ptr, new_offset, factor.len, factor.sp, factor.ep, factor.sp1, factor.ep1);
        }
    }
    
    /* encode data  */
    coder_pv.encode_block_output(os,bfd,block_id,"RLZP2C-PV");
    coder_zzz.encode_block_output(os,bfd_zzz,block_id,"RLZP2C-ZZZ");
}

template<class t_itr,class t_coder_pv,class t_coder_zzz,class t_stream,class t_idx>
void factor_block_rlzp3(
    size_t block_id,
    t_itr itr,
    size_t block_size,
    t_coder_pv& coder_pv,
    t_coder_zzz& coder_zzz,
    t_stream& os,
    uint64_t page_size_log2,
    std::vector<uint8_t>& preferred_pages,
    t_idx& idx,
    block_factor_data& bfd,
    block_factor_data& bfd_zzz)
{
    /* reset data */
    os.seek(0);
    bfd.reset();
    bfd_zzz.reset();
    for (size_t k = 0; k < preferred_pages.size(); k++) preferred_pages[k] = 0;
    
    /* factorize */
    auto end = itr + block_size;
    auto factor_itr = idx.template factorize<decltype(itr),false>(itr, end);
    size_t syms_encoded = 0;
    while (!factor_itr.finished()) {
        if (factor_itr.len == 0) {
            bfd.add_factor(coder_pv, itr + syms_encoded, 0, 1, 0, 0, 0,0);
            bfd_zzz.add_factor(coder_zzz, itr + syms_encoded, 0, 1, 0, 0, 0,0);
            syms_encoded++;
        } else {
            size_t m = factor_itr.ep - factor_itr.sp + 1;
            /* find preferred in <sp,ep> */
            bool preferred_found = false;
            uint32_t preferred_offset = 0;
            uint32_t smallest_offset = std::numeric_limits<uint32_t>::max();
            for (size_t i = 0; i < m; i++) {
                auto cur_offset = idx.translate_offset(factor_itr.sp + i, factor_itr.len);
                auto cur_page = cur_offset >> page_size_log2;
                if (preferred_pages[cur_page] == 1) {
                    preferred_found = true;
                    preferred_offset = cur_offset;
                    break;
                }
                if(smallest_offset > cur_offset) {
                    smallest_offset = cur_offset;
                }
            }
            /* search <sp1,ep1> if no preferred found */
            bool preferred_found_minus_one = false;
            if(!preferred_found && factor_itr.len > coder_pv.literal_threshold + 1) {
                 size_t m1 = factor_itr.ep1 - factor_itr.sp1 + 1;
                 if(m1 != m) {
                    for (size_t i = 0; i < m1; i++) {
                        auto cur_offset = idx.translate_offset(factor_itr.sp1 + i, factor_itr.len-1);
                        auto cur_page = cur_offset >> page_size_log2;
                        if (preferred_pages[cur_page] == 1) {
                            preferred_found_minus_one = true;
                            preferred_offset = cur_offset;
                            factor_itr.rewind(); // have to go back one symbol in the factorization
                            break;
                        }
                    }
                 }
            }
            
            if(!preferred_found && !preferred_found_minus_one) {
                auto new_page = smallest_offset >> page_size_log2;
                preferred_pages[new_page] = 1;
                preferred_offset = smallest_offset;
            }
            
            auto offset = preferred_offset;
            size_t flen = factor_itr.len;
            if(preferred_found_minus_one) {
                flen--;
            }
            bfd.add_factor(coder_pv, itr + syms_encoded, offset, flen, 0, 0, 0,0);
            bfd_zzz.add_factor(coder_zzz, itr + syms_encoded, offset, flen, 0, 0, 0,0);
            syms_encoded += flen;
        }
        ++factor_itr;
    }
        
    
    /* encode data  */
    coder_pv.encode_block_output(os,bfd,block_id,"RLZP3-PV");
    coder_zzz.encode_block_output(os,bfd_zzz,block_id,"RLZP3-ZZZ");
}


int main(int argc, const char* argv[])
{
    setup_logger(argc, argv);

    /* parse command line */
    LOG(INFO) << "Parsing command line arguments";
    auto args = utils::parse_args(argc, argv);

    /* parse the collection */
    LOG(INFO) << "Parsing collection directory " << args.collection_dir;
    collection col(args.collection_dir);

    /* define RLZ parameters */
    const uint32_t factorization_blocksize = 64 * 1024;
    const uint64_t dict_size_bytes = 256 * 1024 * 1024;
    const uint64_t literal_threshold = 3;
    const uint64_t dict_segment_size_bytes = 1024;
    const uint64_t dict_pointer_width = utils::CLog2<dict_size_bytes>();
    const uint64_t segment_size_log2 = utils::CLog2<dict_segment_size_bytes>();
    const uint64_t dict_page_size = 16 * dict_segment_size_bytes;
    const uint64_t page_size_log2 = utils::CLog2<dict_page_size>();
    const uint64_t in_page_offset_width = utils::CLog2<dict_page_size>();
    const uint64_t page_ptr_width = dict_pointer_width - in_page_offset_width;
    const uint64_t num_pages_in_dict = dict_size_bytes / dict_page_size;

    /* output parameters */
    LOG(INFO) << "factorization_blocksize = " << factorization_blocksize;
    LOG(INFO) << "dict_size_bytes = " << dict_size_bytes;
    LOG(INFO) << "literal_threshold = " << literal_threshold;
    LOG(INFO) << "dict_segment_size_bytes = " << dict_segment_size_bytes;
    LOG(INFO) << "log2(dict_segment_size_bytes) = " << segment_size_log2;
    LOG(INFO) << "dict_page_size = " << dict_page_size;
    LOG(INFO) << "log2(dict_page_size) = " << page_size_log2;
    LOG(INFO) << "num_pages_in_dict = " << num_pages_in_dict;
    LOG(INFO) << "dict_pointer_width = " << dict_pointer_width;
    LOG(INFO) << "in_page_offset_width = " << in_page_offset_width;
    LOG(INFO) << "page_ptr_width = " << page_ptr_width;


    // /* define RLZ types */
    const bool default_search_local_context = false;
    using dict_creation_strategy = dict_uniform_sample_budget_sep<dict_segment_size_bytes>;
    using dict_pruning_strategy = dict_prune_none;
    using factor_selection_strategy = factor_select_first;
    using block_map_type = block_map_uncompressed;
    using csa_type_ = sdsl::csa_wt<sdsl::wt_huff<sdsl::bit_vector_il<64> >, 1, 4096>;
    using dict_index_type = dict_index_csa<csa_type_>;

    using baseline_factor_coder_packed = factor_coder_blocked<literal_threshold, coder::fixed<8>, coder::fixed<dict_pointer_width>, coder::vbyte>;
    using baseline_factor_coder_zlib = factor_coder_blocked<literal_threshold, coder::zlib<9>, coder::zlib<9>, coder::zlib<9> >;

    // using baseline_type_packed = rlz_store_static<dict_creation_strategy,
    //                                               dict_pruning_strategy,
    //                                               dict_index_type,
    //                                               factorization_blocksize,
    //                                               default_search_local_context,
    //                                               factor_selection_strategy,
    //                                               baseline_factor_coder_packed,
    //                                               block_map_type>; 

    using baseline_type_zlib = rlz_store_static<dict_creation_strategy,
                                                  dict_pruning_strategy,
                                                  dict_index_type,
                                                  factorization_blocksize,
                                                  default_search_local_context,
                                                  factor_selection_strategy,
                                                  baseline_factor_coder_zlib,
                                                  block_map_type>; 

    using www_type_zlib = rlz_store_static<dict_local_coverage_norms<dict_segment_size_bytes>,
                                                  dict_pruning_strategy,
                                                  dict_index_type,
                                                  factorization_blocksize,
                                                  default_search_local_context,
                                                  factor_selection_strategy,
                                                  baseline_factor_coder_zlib,
                                                  block_map_type>; 

    using p0pv_coder_t = factor_coder_blocked_subdict_pv<literal_threshold,
                                                  dict_page_size,
                                                  num_pages_in_dict,
                                                  dict_size_bytes,
                                                  coder::elias_fano>;    

    using p0zzz_coder_t = factor_coder_blocked_subdict_zzz<literal_threshold,
                                                  dict_page_size,
                                                  num_pages_in_dict,
                                                  dict_size_bytes,
                                                  coder::elias_fano>;         
                                                  
    using p1pv_coder_t = factor_coder_blocked_subdict_pv<literal_threshold,
                                                  dict_page_size,
                                                  num_pages_in_dict,
                                                  dict_size_bytes,
                                                  coder::interpolative>;    

    using p1zzz_coder_t = factor_coder_blocked_subdict_zzz<literal_threshold,
                                                  dict_page_size,
                                                  num_pages_in_dict,
                                                  dict_size_bytes ,
                                                  coder::interpolative>;  
                                                  
    using factor_coder_zlib_transposed = factor_coder_blocked<literal_threshold, coder::zlib<9>, coder::zlibt<9>, coder::zlib<9> >;

    {
        // auto rlz_store = baseline_type_zlib::builder{}
        //                      .set_rebuild(args.rebuild)
        //                      .set_threads(args.threads)
        //                      .set_dict_size(dict_size_bytes)
        //                      .build_or_load(col);
     
        auto rlz_store = www_type_zlib::builder{}
                             .set_rebuild(args.rebuild)
                             .set_threads(args.threads)
                             .set_dict_size(dict_size_bytes)
                             .build_or_load(col);
     
        LOG(INFO) << "Load dict index";
        dict_index_type dict_idx(col, false);
     
        LOG(INFO) << "Create or load segment table";
        auto segment_table = create_or_load_segment_table(col,rlz_store,dict_idx,dict_segment_size_bytes,dict_size_bytes);
	     
        
        LOG(INFO) << "name;"
                  << "block_id;"
                  << "len_size_bytes;"
                  << "literal_size_bytes;"
                  << "offset_size_bytes;"
                  << "page_table_size_bytes;"
                  << "num_factors;"
                  << "num_literals;"
                  << "num_offsets;"
                  << "num_pages_in_block;"
                  << "use_page_table;"
                  << "total_size_bytes";

        std::vector<factor_data> factor_data;
        block_factor_data bfd(factorization_blocksize);
        block_factor_data bfd_zzz(factorization_blocksize);
        sdsl::bit_vector tmp_buffer(1024*1024*8*32);
        bit_ostream<sdsl::bit_vector> fs(tmp_buffer);
        p0pv_coder_t p0pv_coder; p0zzz_coder_t p0zzz_coder;
        p1pv_coder_t p1pv_coder; p1zzz_coder_t p1zzz_coder;
        baseline_factor_coder_packed pv_coder; baseline_factor_coder_zlib zzz_coder;
        factor_coder_zlib_transposed zzzt_coder;
        std::vector<uint8_t> preferred_pages(num_pages_in_dict+1);

        // auto fitr = rlz_store.factors_begin();
        // auto fend = rlz_store.factors_end();
        // while(fitr != fend) {
        //     const auto& f = *fitr;
        //     factor_data.push_back(f);
        //     if (fitr.last_factor_in_block()) {
        //         // compress_rlz_block_inmem("RLZ-PV",pv_coder,fitr.block_id,fs,factor_data,bfd);
        //         // compress_rlz_block_inmem("RLZ-ZZZ",zzz_coder,fitr.block_id,fs,factor_data,bfd);
        //         // compress_rlz_block_inmem("RLZ-ZZZT",zzzt_coder,fitr.block_id,fs,factor_data,bfd);
        //         // compress_rlz_block_inmem("RLZP0-PV",p0pv_coder,fitr.block_id,fs,factor_data,bfd);
        //         // compress_rlz_block_inmem("RLZP0-ZZZ",p0zzz_coder,fitr.block_id,fs,factor_data,bfd);
        //         // compress_rlz_block_inmem("RLZP1-PV",p1pv_coder,fitr.block_id,fs,factor_data,bfd);
        //         // compress_rlz_block_inmem("RLZP1-ZZZ",p1zzz_coder,fitr.block_id,fs,factor_data,bfd);
        //         // compress_rlzp2a_block_inmem(p1pv_coder,p1zzz_coder,
        //         //     fitr.block_id,fs,factor_data,segment_table,segment_size_log2,
        //         // //     dict_idx,bfd,bfd_zzz);
        //         // compress_rlzp2b_block_inmem(p1pv_coder,p1zzz_coder,
        //         //     fitr.block_id,fs,factor_data,segment_table,segment_size_log2,
        //         //     page_size_log2,preferred_pages,dict_idx,bfd,bfd_zzz);
        //         // compress_rlzp2c_block_inmem(p1pv_coder,p1zzz_coder,
        //         //     fitr.block_id,fs,factor_data,segment_table,segment_size_log2,
        //         //     page_size_log2,preferred_pages,dict_idx,bfd,bfd_zzz);
        //         // compress_rlzp3_block_inmem(p1pv_coder,p1zzz_coder,
        //         //     fitr.block_id,fs,factor_data,segment_table,segment_size_log2,page_size_log2,preferred_pages,dict_idx,bfd);
        //         factor_data.clear();
        //     }
        //     ++fitr;
        // }
        
        // /* P3 */
        LOG(INFO) << "P3 Write dict";
        auto new_dict_file_name = dict_segment_score<dict_segment_size_bytes,dict_page_size>::file_name(col, dict_size_bytes);
        if (!utils::file_exists(new_dict_file_name)) {
            /* (1) get the old dict */
            const auto& existing_dict = rlz_store.dict;
            auto existing_dict_hash = rlz_store.m_dict_hash;
            col.file_map[KEY_DICT] = new_dict_file_name;
            auto new_dict = sdsl::write_out_buffer<8>::create(col.file_map[KEY_DICT]);
            /* (4) write the reordered dictionary */
            for (size_t i = 0; i < segment_table.size(); i++) {
                auto old_segment_id = segment_table[i];
                auto old_segment_itr = existing_dict.begin() + (old_segment_id * dict_segment_size_bytes);
                for (size_t j = 0; j < dict_segment_size_bytes; j++) {
                    new_dict.push_back(*old_segment_itr++);
                }
            }
            new_dict.push_back(0); // EOD symbol
        }
        else {
            col.file_map[KEY_DICT] = new_dict_file_name;
        }
        col.compute_dict_hash();
        LOG(INFO) << "Create dict index over segment score sorted dict";
        dict_index_type ss_dict_idx(col, false);
         
        {
            
            const sdsl::int_vector_mapped_buffer<8> text(col.file_map[KEY_TEXT]);
            size_t n = text.size();
            size_t num_blocks = n / factorization_blocksize;
            auto itr = text.begin();
            for(size_t i=0;i<num_blocks;i++) {
                factor_block_rlzp3(i,itr,factorization_blocksize,p1pv_coder,p1zzz_coder,
                    fs,page_size_log2,preferred_pages,ss_dict_idx,bfd,bfd_zzz);
                itr += factorization_blocksize;
            }
            itr = itr + num_blocks*factorization_blocksize;
            auto left = n % factorization_blocksize;
            if(left) {
                factor_block_rlzp3(num_blocks,itr,left,p1pv_coder,p1zzz_coder,
                    fs,page_size_log2,preferred_pages,ss_dict_idx,bfd,bfd_zzz);
            }
        }
     
    }
    
    return EXIT_SUCCESS;
}
