
#include "count_min_sketch.hpp"
#include "sdsl/int_vector.hpp"

#include <cstdio>
#include <fstream>

int main(int argc,const char* argv[])
{
    if(argc!=2) return -1;
    auto file_name = argv[1];
    FILE* fp = std::fopen(file_name, "r");
    uint64_t cnts[256*256] = {0};
    size_t k=10;
    count_min_sketch_topk<uint64_t> cms(k);
    int cc = std::getc(fp);
    int c;
    auto start = std::chrono::high_resolution_clock::now();
    while ((c = std::getc(fp)) != EOF) {
        uint64_t sym = (cc << 8) + c;
        cms.update(sym);
        cnts[sym]++;
        cc = c;
    }
    auto stop = std::chrono::high_resolution_clock::now();
    std::fclose(fp);

    std::cout << "epsilon = " << cms.sketch.epsilon << std::endl;
    std::cout << "delta   = " << cms.sketch.delta << std::endl;
    std::cout << "w       = " << cms.sketch.w << std::endl;
    std::cout << "d       = " << cms.sketch.d << std::endl;
    std::cout << "size in bytes = " << cms.sketch.size_in_bytes() << std::endl;
    std::cout << "estimation_error = " << cms.sketch.estimation_error() << std::endl;
    std::cout << "estimation_prob = " << cms.sketch.estimation_probability() << std::endl;

    // for(size_t i=0;i<256*256;i++) {
    //     auto est = cms.estimate(i);
    //     auto real = cnts[i];
    //     auto diff = est-real;
    //     std::cout << "sym=" << i 
    //               << " real = " << real 
    //               << " est= " << est 
    //               << " diff=" << diff 
    //               << std::endl;
    // }
    std::vector<std::pair<uint64_t,uint64_t>> sorted_cnts;
    for(size_t i=0;i<256*256;i++) {
        sorted_cnts.emplace_back(cnts[i],i);
    }
    std::sort(sorted_cnts.begin(),sorted_cnts.end(),std::greater<std::pair<uint64_t,uint64_t>>());
    for(size_t i=0;i<k;i++) {
        std::cout << "("<<i<<") item = " << sorted_cnts[i].second << " freq = " << sorted_cnts[i].first << std::endl;
    }

    auto topk = cms.topk();
    size_t j=1;
    for(const auto& item : topk) {
        std::cout << "("<<j<<") item = " << item.second << " estimate=" << item.first << " real=" << cnts[item.second] << std::endl;
        j++;
    }

    auto time_spent = stop-start;
    std::cout << "time in sec = " << std::chrono::duration_cast<std::chrono::milliseconds>(time_spent).count() / 1000.0f << std::endl;

    return 0;
}
