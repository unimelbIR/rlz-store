#define ELPP_THREAD_SAFE

#include "utils.hpp"
#include "collection.hpp"
#include "rlz_utils.hpp"

#include "indexes.hpp"

#include "logging.hpp"
INITIALIZE_EASYLOGGINGPP

int main(int argc, const char* argv[])
{
    setup_logger(argc, argv);

    /* parse command line */
    LOG(INFO) << "Parsing command line arguments";
    auto args = utils::parse_args(argc, argv);

    /* parse the collection */
    LOG(INFO) << "Parsing collection directory " << args.collection_dir;
    collection col(args.collection_dir);

    sdsl::int_vector<8> text;
    sdsl::load_from_file(text,col.file_map[KEY_TEXT]);
    for(size_t i=0;i<text.size();i++) {
        if(text[i]==1) text[i] = 0xFF;        
    }
    sdsl::store_to_file(text,col.file_map[KEY_TEXT]);

    return EXIT_SUCCESS;
}
