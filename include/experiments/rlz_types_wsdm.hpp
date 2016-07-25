#pragma once

#include "indexes.hpp"

#include "wt_flat.hpp"

using www_csa_type = sdsl::csa_wt<wt_flat<sdsl::bit_vector>, 1, 4096>;
using csa_type = sdsl::csa_wt<sdsl::wt_huff<sdsl::bit_vector_il<64> >, 4, 4096>;
const uint32_t www_uniform_sample_block_size = 1024;
const uint32_t factorization_blocksize = 1024*1024;


using rlz_store_static_single = rlz_store_static<dict_local_coverage_norms<1024,16,512,std::ratio<1,2>>,
                                dict_prune_none,
                                dict_index_csa<www_csa_type>, //to update to sa
                                factorization_blocksize,
                                false,
                                factor_select_first,
                                factor_coder_blocked<3, coder::zlib<9>, coder::zlib<9>, coder::zlib<9>>,
                                block_map_uncompressed>;

using rlz_store_static_multi =  rlz_store_static<dict_multibale_local_coverage_norms<1024,16,512,std::ratio<1,2>>,
                                dict_prune_none,
                                dict_index_csa<www_csa_type>, //to update to sa
                                factorization_blocksize,
                                false,
                                factor_select_first,
                                factor_coder_blocked<3, coder::zlib<9>, coder::zlib<9>, coder::zlib<9>>,
                                block_map_uncompressed>;

