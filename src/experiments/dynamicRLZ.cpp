// syntax dynamicRLZ #collectionPath(-c) #bales(-b) #window_size(-w) #dic_size(-d) #mode(-m)

#define ELPP_THREAD_SAFE
#define ELPP_STL_LOGGING

#include "utils.hpp"
#include "collection.hpp"
#include "rlz_utils.hpp"
#include "count_min_sketch.hpp"
#include "chunk_freq_estimator.hpp"
#include <cstdio>
#include <unordered_set>
#include "logging.hpp"
#include "wt_flat.hpp"
#include <algorithm>
#include <iostream>
#include <fstream>
#include "indexes.hpp"

//to change to faster factorization by Matt
using www_csa_type = sdsl::csa_wt<wt_flat<sdsl::bit_vector>, 1, 4096>;
using csa_type = sdsl::csa_wt<sdsl::wt_huff<sdsl::bit_vector_il<64> >, 4, 4096>;

INITIALIZE_EASYLOGGINGPP

//function to load k-mers
void addHistMers(std::unordered_set<uint64_t> &history_mers, std::string inFile) {
    auto start = hrclock::now();

    //read into data structure 
    sdsl::read_only_mapper<8> dict(inFile); 
    fixed_hasher<16> rk;

    for(size_t j=0;j<dict.size();j++) {
        auto sym = dict[j];
        auto hash = rk.update(sym);
        if(j < 16-1) continue; 
        else {
            history_mers.emplace(hash);
        }
    }
    LOG(INFO) << "\t" << "history_mers size: " << history_mers.size();
	 
    // //write to output file
    // LOG(INFO) << "\t" << "Writing out history mers...";
    // auto history_out = sdsl::write_out_buffer<64>::create(file_history_mers);
    // {
    //     auto itr = history_mers.begin();
    //     while(itr != history_mers.end()) {
    //        //std::copy(itr, ++itr, std::back_inserter(history_out));
    //        // itr++;
    // 	   history_out.push_back(*itr);
    // 	   itr++;
    //     }
    // }
    auto stop = hrclock::now();    
    LOG(INFO) << "\t" << "History  time = " << duration_cast<milliseconds>(stop-start).count() / 1000.0f << " sec";
}

void combineDicts(const std::vector<std::string>& file_dicts, std::string file_combined, int start) {
    //read and write
    auto begin = hrclock::now();
    auto out = sdsl::write_out_buffer<8>::create(file_combined);

    //read into data structure
    for(size_t i=start;i<file_dicts.size();i++) {
        sdsl::read_only_mapper<8> dict(file_dicts[i]); 

        for(size_t j=0;j<dict.size()-1;j++) {
            out.push_back(dict[j]);
        }
    }
    out.push_back('\0');
    auto stop = hrclock::now();    
    LOG(INFO) << "\t" << "Combining dictionary time = " << duration_cast<milliseconds>(stop-begin).count() / 1000.0f << " sec";
}

template<class t_idx>
size_t zlib_dict_bits (t_idx& idx) {
    size_t bits_compressed_dict = 0;
    {
        const uint8_t* dict = (const uint8_t*) idx.dict.data();
        size_t dict_size = idx.dict.size();
        std::vector<uint8_t> dict_buf(dict_size*2);
        uint8_t* out_buf = dict_buf.data();
        size_t out_len = dict_buf.size();
        int cok = compress2(out_buf,&out_len,dict,dict_size,9);
        if(cok != Z_OK) {
            LOG(FATAL) << "error compressing dictionary.";
        }
        bits_compressed_dict = out_len * 8;
    }
    return bits_compressed_dict;
}

template<class t_idx>
void compute_archive_ratio(collection& col, std::ofstream &out, t_idx& idx, size_t bits_compressed_dict, std::string index_name)
{
    size_t bits_encoding = idx.factor_text.size();
    size_t bits_blockmap = idx.block_map.size_in_bytes() * 8;
    size_t encoding_bits_total = bits_encoding + bits_compressed_dict + bits_blockmap;
    size_t text_size_bits = idx.size() * 8;
    
    double archive_ratio = 100.0 * double(encoding_bits_total) / double(text_size_bits);
    
    //output to a separate file
    out << index_name << " Encoding size = " << double(bits_encoding) / double(1024*1024*8) << " MiB" << std::endl;
    out << index_name << " Blockmap size = " << double(bits_blockmap) / double(1024*1024*8) << " MiB" << std::endl;
    out << index_name << " Compressed dict size = " << double(bits_compressed_dict) / double(1024*1024*8) << " MiB"  << std::endl;
    out << index_name << " Archive Ratio = " << archive_ratio << " %"  << std::endl << std::endl;

    LOG(INFO) << index_name << " Encoding size = " << double(bits_encoding) / double(1024*1024*8) << " MiB";
    LOG(INFO) << index_name << " Blockmap size = " << double(bits_blockmap) / double(1024*1024*8) << " MiB";
    LOG(INFO) << index_name << " Compressed dict size = " << double(bits_compressed_dict) / double(1024*1024*8) << " MiB";
    LOG(INFO) << index_name << " Archive Ratio = " << archive_ratio << " %";
}

//return the compressed dict size
size_t create_indexes_combine(collection& col, size_t dict_size_in_bytes, int ctype, std::ofstream &out, utils::cmdargs_t& args, bool isCombinedDict, size_t bits_compressed_dict = 0)
{    /* create rlz index */
    const uint32_t factorization_blocksize = 64 * 1024;
    {
        auto rlz_store = rlz_store_static<dict_local_coverage_norms<1024,16,512,std::ratio<1,2>>,
                    dict_prune_none,
                    dict_index_csa<www_csa_type>, //to update to sa
                    factorization_blocksize,
                    false,
                    factor_select_first,
                    factor_coder_blocked<3, coder::zlib<9>, coder::zlib<9>, coder::zlib<9>>,
                    block_map_uncompressed>::builder{}
                    .set_rebuild(args.rebuild)
                    .set_threads(args.threads)
                    .set_dict_size(dict_size_in_bytes)
                    .build_or_load(col, NULL, ctype);
                    
        uint32_t text_size_mib = rlz_store.size() / (1024*1024);
        uint32_t dict_size_mib = dict_size_in_bytes / (1024*1024);
        std::string index_name = "GOV2S-WWW-COMBINE-" + std::to_string(text_size_mib) + "-C" + col.text + "-rw" + std::to_string(ctype) + "-" + std::to_string(dict_size_mib);
        verify_index(col, rlz_store);
        size_t store_bits_compressed_dict = zlib_dict_bits(rlz_store);
        if(isCombinedDict)
            compute_archive_ratio(col,out,rlz_store,bits_compressed_dict,index_name); //not using compressed combined dict size
        else
            compute_archive_ratio(col,out,rlz_store,store_bits_compressed_dict,index_name);
        return store_bits_compressed_dict;
    }
}

//return the compressed dict size
size_t create_indexes_cascade(collection& col, size_t dict_size_in_bytes, int ctype, std::ofstream &out, std::unordered_set<uint64_t> *history_mers, utils::cmdargs_t& args, bool toFactorize)
{    /* create rlz index */
    const uint32_t factorization_blocksize = 64 * 1024;
    {
        auto rlz_store = rlz_store_static<dict_multibale_local_coverage_norms<1024,16,512,std::ratio<1,2>>,
                    dict_prune_none,
                    dict_index_csa<www_csa_type>, //to update to sa
                    factorization_blocksize,
                    false,
                    factor_select_first,
                    factor_coder_blocked<3, coder::zlib<9>, coder::zlib<9>, coder::zlib<9>>,
                    block_map_uncompressed>::builder{}
                    .set_rebuild(args.rebuild)
                    .set_threads(args.threads)
                    .set_dict_size(dict_size_in_bytes)
                    .build_or_load(col, history_mers, ctype, toFactorize);

        size_t store_bits_compressed_dict = zlib_dict_bits(rlz_store);

        if(toFactorize) { //if factorization needed
            uint32_t text_size_mib = rlz_store.size() / (1024*1024);
            uint32_t dict_size_mib = dict_size_in_bytes / (1024*1024);
            std::string index_name = "GOV2S-WWW-CASCADE-"  + std::to_string(text_size_mib) + "-C" + col.text + "-rw" + std::to_string(ctype) + "-" + std::to_string(dict_size_mib);
            verify_index(col, rlz_store);
            // if(isCombinedDict)
            //     compute_archive_ratio(col,out,rlz_store,bits_compressed_dict,index_name); //not using compressed combined dict size
            // else
                compute_archive_ratio(col,out,rlz_store,store_bits_compressed_dict,index_name);
        }
        return store_bits_compressed_dict;
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
    auto n = args.num_bales;
    auto w = args.window_size;
    auto b = args.test_bale; // b < n
    auto mode = args.mode;
    auto dict_size = args.dict_size_in_bytes;
    auto rebuild = args.rebuild;
    std::vector<std::string> dicts; //store dict file names
    std::string out_file; //store the final combined dict file name
    if(b >= n) {
        LOG(INFO) << "\t" << "TEST BALE does not exist! ";
        return -1; //error and terminate
    }
    //setup separate results output
    std::ofstream out(mode + "-b" + std::to_string(b) + "-w" + std::to_string(w) + "-d" + std::to_string(dict_size/(1024*1024)) + ".out", std::ios::app);
    if(!out) {
        std::cout << "error opening file" << std::endl;
        return -1;
    }
    LOG(INFO) << "\t" << "Entering Dynamic Mode = " << mode;

    //simple combine mode
    if(mode == "combine") {
        size_t combined_dict_size_compressed = 0; //will store the test bale dict size in bits
        for (int i = 0; i <= b; i++) {
            collection col(args.collection_dir, std::to_string(i));
            combined_dict_size_compressed = create_indexes_combine(col,dict_size,0,out,args,false);
            dicts.push_back(col.file_map[KEY_DICT]);
        }
        //combine setup     
        auto start = std::max(0,b-w);
        auto real_w = b-start;
        auto c_size = dict_size * (real_w + 1);
        collection col(args.collection_dir, std::to_string(b));
        std::string out_file = dict_local_coverage_norms<1024,16,512,std::ratio<1,2>>::dict_file_name(col, c_size, real_w);

        out << "Finally......" << std::endl;
        out << "Combining simple dictionaries for Bale = " << b << std::endl;
        out << "Dictionary Size in use = " << std::to_string(dict_size/(1024*1024)) << std::endl;
        out << "Context Size = " << w << std::endl;
        out << "Real test bale Context Size = " << real_w << std::endl;

        LOG(INFO) << "\t" << "Finally......";
        LOG(INFO) << "\t" << "Combining simple dictionaries for Bale = " << b;
        LOG(INFO) << "\t" << "Dictionary Size in use = " << std::to_string(dict_size/(1024*1024));
        LOG(INFO) << "\t" << "Context Size = " << w;
        LOG(INFO) << "\t" << "Real Context Size = " << real_w;

        if(! utils::file_exists(out_file) || rebuild )
            combineDicts(dicts, out_file, start);

        //factorize for results
        create_indexes_combine(col,c_size,real_w,out,args,true, combined_dict_size_compressed);  //factorise for compression results
    } 

    //more complicated cascade mode
    if(mode == "cascade") {    
        size_t combined_dict_size_compressed = 0; //will store the test bale dict size in bits
        for (int i = 0; i <= b; i++) {
            auto hist_start = std::max(0,i-w);
            std::unordered_set<uint64_t> history_mers;
            history_mers.max_load_factor(0.1); //make faster by losing memory

            //preload history mers if there is any
            for (int h = hist_start; h <= i - 1; h++) {
                //history dicts should already exist and mers should be preloaded
                auto start = std::max(0,h-w);
                auto real_w = h-start;

                //get history dict name
                collection col(args.collection_dir, std::to_string(h));
                std::string c_type = "-rw"+ std::to_string(real_w);

                std::string hist_file = dict_multibale_local_coverage_norms<1024,16,512,std::ratio<1,2>>::dict_file_name(col, dict_size, real_w);
                if(utils::file_exists(hist_file))
                    addHistMers(history_mers, hist_file);
                else {
                    LOG(INFO) << "\t" << "History File: " << hist_file << "do not exist! Program exit!";
                    return -1; //exit(-1);?
                }
            }
            
            //build own dict
            collection col(args.collection_dir, std::to_string(i));
            auto real_w = i-hist_start;
            // std::string c_type = "-rw"+ std::to_string(real_w);
            //don't wanna factorise now!
            if(history_mers.empty())
                combined_dict_size_compressed = create_indexes_cascade(col,dict_size,real_w,out,NULL,args,false);
            else
                combined_dict_size_compressed = create_indexes_cascade(col,dict_size,real_w,out,&history_mers,args,false);
            
            dicts.push_back(col.file_map[KEY_DICT]);
        }

        //combine setup     
        auto start = std::max(0,b-w);
        auto real_w = b-start;
        auto c_size = dict_size * (real_w + 1);
        collection col(args.collection_dir, std::to_string(b));
        std::string c_type = "-rw"+ std::to_string(real_w);
        std::string out_file = dict_multibale_local_coverage_norms<1024,16,512,std::ratio<1,2>>::dict_file_name(col, c_size, real_w);

        out << "Finally......" << std::endl;
        out << "Combining cascaded dictionaries for Bale = " << b << std::endl;
        out << "Dictionary Size in use = " << std::to_string(dict_size/(1024*1024)) << std::endl;
        out << "Context Size = " << w << std::endl;
        out << "Real test bale Context Size = " << real_w << std::endl;

        LOG(INFO) << "\t" << "Finally......";
        LOG(INFO) << "\t" << "Combining cascaded dictionaries for Bale = " << b;
        LOG(INFO) << "\t" << "Dictionary Size in use = " << std::to_string(dict_size/(1024*1024));
        LOG(INFO) << "\t" << "Context Size = " << w;
        LOG(INFO) << "\t" << "Real Context Size = " << real_w;

        if(! utils::file_exists(out_file) || rebuild )
            combineDicts(dicts, out_file, start);

        //factorize for results: warning: may overwrite exsiting combined dicts!
        create_indexes_combine(col,c_size,real_w,out,args,true, combined_dict_size_compressed); //factorise for compression results
    }

    /* create rlz indices */
    // create_indexes<1*1024*1024>(col,args);
    // create_indexes<2*1024*1024>(col,args);
    // create_indexes<4*1024*1024>(col,args);
    // create_indexes<8*1024*1024>(col,args);
    // create_indexes<16*1024*1024>(col,args);
    // create_indexes<32*1024*1024>(col,args);
    // create_indexes<64*1024*1024>(col,args);
    return EXIT_SUCCESS;
}
    
