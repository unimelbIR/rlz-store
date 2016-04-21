#define ELPP_THREAD_SAFE
#define ELPP_STL_LOGGING

#include "utils.hpp"
#include "collection.hpp"
#include "rlz_utils.hpp"

#include "indexes.hpp"

#include "logging.hpp"

INITIALIZE_EASYLOGGINGPP

//generate dict history mers for prefilling
//format: argv[argc-1] is the output histroy mers, other inputs by input dict order
int main(int argc, const char* argv[])
{
    //assign string labels
    auto file_combined = argv[argc-1];  
    std::string file_dicts[argc-2];

    for(int i=1;i<argc-1;i++) { 
        file_dicts[i-1] = argv[i];
    }        

    //read and write
    auto start = hrclock::now();
    auto out = sdsl::write_out_buffer<8>::create(file_combined);

    //read into data structure
    LOG(INFO) << "\t" << "Reading dictionary files...";
    for(int i=0;i<argc-2;i++) {
        sdsl::read_only_mapper<8> dict(file_dicts[i]); 

        for(size_t j=0;j<dict.size()-1;j++) {
            out.push_back(dict[j]);
        }
    }
    out.push_back('\0');
    auto stop = hrclock::now();    
    LOG(INFO) << "\t" << "Combining dictionary time = " << duration_cast<milliseconds>(stop-start).count() / 1000.0f << " sec";
}
