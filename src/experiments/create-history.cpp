#define ELPP_THREAD_SAFE
#define ELPP_STL_LOGGING

#include "utils.hpp"
#include "collection.hpp"
#include "rlz_utils.hpp"

#include "indexes.hpp"

#include "logging.hpp"

#include "count_min_sketch.hpp"
#include "chunk_freq_estimator.hpp"

#include <unordered_set>

INITIALIZE_EASYLOGGINGPP

//generate dict history mers for prefilling
//format: argv[argc-1] is the output histroy mers, other inputs by input dict order
int main(int argc, const char* argv[])
{
    //assign string labels
    auto file_history_mers = argv[argc-1];  
    std::string file_dicts[argc-2];
    const uint32_t mers_size = 16; //hard code

    for(int i=1;i<argc-1;i++) { 
        file_dicts[i-1] = argv[i];
    }        

    //read and write history mers
    std::unordered_set<uint64_t> history_mers;
    history_mers.max_load_factor(0.1); //make faster by losing memory
    auto start = hrclock::now();

    //read into data structure
    LOG(INFO) << "\t" << "Reading dictionary files...";
    for(int i=0;i<argc-2;i++) {
        sdsl::read_only_mapper<8> dict(file_dicts[i]); 
        fixed_hasher<mers_size> rk;

        for(size_t j=0;j<dict.size();j++) {
            auto sym = dict[j];
            auto hash = rk.update(sym);
            if(j < mers_size-1) continue; 
            else {
                history_mers.emplace(hash);
            }
        }
    }
     LOG(INFO) << "\t" << "history_mers size: " << history_mers.size();
	 
    //write to output file
    LOG(INFO) << "\t" << "Writing out history mers...";
    auto history_out = sdsl::write_out_buffer<64>::create(file_history_mers);
    {
        auto itr = history_mers.begin();
        while(itr != history_mers.end()) {
           //std::copy(itr, ++itr, std::back_inserter(history_out));
           // itr++;
    	   history_out.push_back(*itr);
    	   itr++;
        }
    }
    auto stop = hrclock::now();    
    LOG(INFO) << "\t" << "History  time = " << duration_cast<milliseconds>(stop-start).count() / 1000.0f << " sec";
}
