#define ELPP_THREAD_SAFE
#define ELPP_STL_LOGGING

#include "utils.hpp"
#include "collection.hpp"
#include "rlz_utils.hpp"

#include "indexes.hpp"

#include "logging.hpp"
INITIALIZE_EASYLOGGINGPP


template <class t_idx>
void bench_index_rand(collection& col,const t_idx& idx,size_t dict_size_in_bytes,std::string name)
{
	utils::flush_cache();

	uint64_t blocks_to_decode = 10000ULL;
	uint64_t block_ret_size = 64*1024;

	std::vector<uint64_t> byte_offsets(blocks_to_decode);
	std::mt19937 g(1234);
	std::uniform_int_distribution<uint64_t> dis(0, idx.text_size - 1 - block_ret_size);
	for(size_t i=0;i<blocks_to_decode;i++) {
		byte_offsets[i] = dis(g);
	}
	
	std::vector<uint8_t> block_content(idx.encoding_block_size);
	block_factor_data bfd(idx.encoding_block_size);

	auto itr = idx.begin();
	auto start = hrclock::now();
	size_t checksum = 0;
	size_t num_syms = 0;
	for(const auto bo: byte_offsets) {
		auto text_ret_offset = bo;
		itr.seek(text_ret_offset);
		num_syms += block_ret_size;
		for (size_t i = 0; i < block_ret_size; i++) {
			checksum += *itr;
			++itr;
		}
	}
    auto stop = hrclock::now();
    if (checksum == 0) {
        LOG(ERROR) << "text decoding speed checksum error";
    }


    auto time_ms = duration_cast<milliseconds>(stop - start);
    LOG(INFO) << name << ";"
    		  << dict_size_in_bytes << ";"
    		  << idx.encoding_block_size << ";"
    		  << time_ms.count() << ";"
    		  << num_syms << ";"
    		  << blocks_to_decode << ";"
    		  << block_ret_size << ";"
    		  << checksum << ";"
              << idx.size_in_bytes() << ";"
              << idx.size() << ";"
              << col.path << ";"
    		  << "RAND";
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
        using rlz_type_t = rlz_store_static<dict_creation_strategy,
                                                  dict_pruning_strategy,
                                                  dict_index_type,
                                                  factorization_blocksize,
                                                  default_search_local_context,
                                                  factor_selection_strategy,
                                                  baseline_factor_coder_zlib,
                                                  block_map_type>; 

        auto rlz_store = rlz_type_t::builder{}
                             .set_rebuild(args.rebuild)
                             .set_threads(args.threads)
                             .set_dict_size(dict_size_bytes)
                             .build_or_load(col);
    
        bench_index_rand(col,rlz_store,dict_size_bytes,"RLZ-ZZZ");
    }
    
    {
        using rlz_type_t = rlz_store_static<dict_creation_strategy,
                                                  dict_pruning_strategy,
                                                  dict_index_type,
                                                  factorization_blocksize,
                                                  default_search_local_context,
                                                  factor_selection_strategy,
                                                  baseline_factor_coder_packed,
                                                  block_map_type>; 

        auto rlz_store = rlz_type_t::builder{}
                             .set_rebuild(args.rebuild)
                             .set_threads(args.threads)
                             .set_dict_size(dict_size_bytes)
                             .build_or_load(col);
    
        bench_index_rand(col,rlz_store,dict_size_bytes,"RLZ-PV");
    }
    
    {
        using rlz_type_t = rlz_store_static<dict_creation_strategy,
                                                  dict_pruning_strategy,
                                                  dict_index_type,
                                                  factorization_blocksize,
                                                  default_search_local_context,
                                                  factor_selection_strategy,
                                                  p0pv_coder_t,
                                                  block_map_type>; 

        auto rlz_store = rlz_type_t::builder{}
                             .set_rebuild(args.rebuild)
                             .set_threads(args.threads)
                             .set_dict_size(dict_size_bytes)
                             .build_or_load(col);
    
        bench_index_rand(col,rlz_store,dict_size_bytes,"RLZP0-PV");
    }
    
    {
        using rlz_type_t = rlz_store_static<dict_creation_strategy,
                                                  dict_pruning_strategy,
                                                  dict_index_type,
                                                  factorization_blocksize,
                                                  default_search_local_context,
                                                  factor_selection_strategy,
                                                  p0zzz_coder_t,
                                                  block_map_type>; 

        auto rlz_store = rlz_type_t::builder{}
                             .set_rebuild(args.rebuild)
                             .set_threads(args.threads)
                             .set_dict_size(dict_size_bytes)
                             .build_or_load(col);
    
        bench_index_rand(col,rlz_store,dict_size_bytes,"RLZP0-ZZZ");
    }
    
    {
        using rlz_type_t = rlz_store_static<dict_creation_strategy,
                                                  dict_pruning_strategy,
                                                  dict_index_type,
                                                  factorization_blocksize,
                                                  default_search_local_context,
                                                  factor_selection_strategy,
                                                  p1pv_coder_t,
                                                  block_map_type>; 

        auto rlz_store = rlz_type_t::builder{}
                             .set_rebuild(args.rebuild)
                             .set_threads(args.threads)
                             .set_dict_size(dict_size_bytes)
                             .build_or_load(col);
    
        bench_index_rand(col,rlz_store,dict_size_bytes,"RLZP1-PV");
    }
    
    {
        using rlz_type_t = rlz_store_static<dict_creation_strategy,
                                                  dict_pruning_strategy,
                                                  dict_index_type,
                                                  factorization_blocksize,
                                                  default_search_local_context,
                                                  factor_selection_strategy,
                                                  p1zzz_coder_t,
                                                  block_map_type>; 

        auto rlz_store = rlz_type_t::builder{}
                             .set_rebuild(args.rebuild)
                             .set_threads(args.threads)
                             .set_dict_size(dict_size_bytes)
                             .build_or_load(col);
    
        bench_index_rand(col,rlz_store,dict_size_bytes,"RLZP1-ZZZ");
    }
    
    
    return EXIT_SUCCESS;
}
