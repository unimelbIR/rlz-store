#define ELPP_THREAD_SAFE
#define ELPP_STL_LOGGING

#include "utils.hpp"
#include "collection.hpp"
#include "rlz_utils.hpp"

#include "indexes.hpp"

#include "logging.hpp"
INITIALIZE_EASYLOGGINGPP


template <class t_idx>
void output_offset_stats(t_idx& idx, std::string name)
{
    for (size_t block_id = 0; block_id < idx.block_map.num_blocks(); block_id++) {
        // factor stats
        auto fdat = idx.block_factors(block_id);
        auto block_size_data = fdat.first;
        LOG(INFO) << name << ";"
                  << block_id << ";"
                  << block_size_data.offset_bytes;
    }
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
    const uint64_t dict_page_size = 16 * dict_segment_size_bytes;
    const uint64_t in_page_offset_width = utils::CLog2<dict_page_size>();
    const uint64_t page_ptr_width = dict_pointer_width - in_page_offset_width;
    const uint64_t num_pages_in_dict = dict_size_bytes / dict_page_size;

    /* output parameters */
    LOG(INFO) << "factorization_blocksize = " << factorization_blocksize;
    LOG(INFO) << "dict_size_bytes = " << dict_size_bytes;
    LOG(INFO) << "literal_threshold = " << literal_threshold;
    LOG(INFO) << "dict_segment_size_bytes = " << dict_segment_size_bytes;
    LOG(INFO) << "dict_page_size = " << dict_page_size;
    LOG(INFO) << "num_pages_in_dict = " << num_pages_in_dict;
    LOG(INFO) << "dict_pointer_width = " << dict_pointer_width;
    LOG(INFO) << "in_page_offset_width = " << in_page_offset_width;
    LOG(INFO) << "page_ptr_width = " << page_ptr_width;


    // /* define RLZ types */
    const bool default_search_local_context = false;
    using dict_creation_strategy = dict_uniform_sample_budget<dict_segment_size_bytes>;
    using dict_pruning_strategy = dict_prune_none;
    using factor_selection_strategy = factor_select_first;
    using block_map_type = block_map_uncompressed;
    using csa_type_ = sdsl::csa_wt<sdsl::wt_huff<sdsl::bit_vector_il<64> >, 1, 4096>;
    using dict_index_type = dict_index_csa<csa_type_>;

    using baseline_factor_coder_packed = factor_coder_blocked<literal_threshold, coder::fixed<8>, coder::fixed<dict_pointer_width>, coder::vbyte>;
    using baseline_factor_coder_zlib = factor_coder_blocked<literal_threshold, coder::zlib<9>, coder::zlib<9>, coder::zlib<9> >;

    using baseline_type_packed = rlz_store_static<dict_creation_strategy,
                                                  dict_pruning_strategy,
                                                  dict_index_type,
                                                  factorization_blocksize,
                                                  default_search_local_context,
                                                  factor_selection_strategy,
                                                  baseline_factor_coder_packed,
                                                  block_map_type>; 

    using baseline_type_zlib = rlz_store_static<dict_creation_strategy,
                                                  dict_pruning_strategy,
                                                  dict_index_type,
                                                  factorization_blocksize,
                                                  default_search_local_context,
                                                  factor_selection_strategy,
                                                  baseline_factor_coder_zlib,
                                                  block_map_type>; 

    using p0_coder = factor_coder_blocked_subdict<literal_threshold,dict_page_size,num_pages_in_dict,coder::fixed<8>, coder::elias_fano, coder::fixed<dict_pointer_width>, coder::vbyte >;                                              
    using rlz_p0 = rlz_store_static<dict_creation_strategy,
                                                  dict_pruning_strategy,
                                                  dict_index_type,
                                                  factorization_blocksize,
                                                  default_search_local_context,
                                                  factor_selection_strategy,
                                                  baseline_factor_coder_zlib,
                                                  block_map_type>; 

    {
        auto rlz_store = baseline_type_packed::builder{}
                             .set_rebuild(args.rebuild)
                             .set_threads(args.threads)
                             .set_dict_size(dict_size_bytes)
                             .build_or_load(col);

        verify_index(col, rlz_store);
        output_offset_stats(rlz_store,"PACKED_BINARY_INTEGERS");
    }
    {
        auto rlz_store = baseline_type_zlib::builder{}
                             .set_rebuild(args.rebuild)
                             .set_threads(args.threads)
                             .set_dict_size(dict_size_bytes)
                             .build_or_load(col);

        verify_index(col, rlz_store);
        output_offset_stats(rlz_store,"ZZZ");
    }

    return EXIT_SUCCESS;
}
