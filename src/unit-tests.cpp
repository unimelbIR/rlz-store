#include "gtest/gtest.h"
#include "sdsl/int_vector.hpp"
#include "bit_streams.hpp"
#include "bit_coders.hpp"
#include "factor_data.hpp"
#include "factor_coder.hpp"
#include <functional>
#include <random>

#include "utils.hpp"

#include "logging.hpp"
INITIALIZE_EASYLOGGINGPP

TEST(bit_stream, vbyte)
{
    size_t n = 20;
    std::mt19937 gen(4711);
    std::uniform_int_distribution<uint64_t> dis(1, 100000);

    for (size_t i = 0; i < n; i++) {
        size_t len = dis(gen);
        std::vector<uint32_t> A(len);
        for (size_t j = 0; j < len; j++)
            A[j] = dis(gen);
        std::sort(A.begin(), A.end());
        auto last = std::unique(A.begin(), A.end());
        auto n = std::distance(A.begin(), last);
        coder::vbyte c;
        sdsl::bit_vector bv;
        {
            bit_ostream<sdsl::bit_vector> os(bv);
            c.encode(os, A.data(), n);
        }
        std::vector<uint32_t> B(n);
        {
            bit_istream<sdsl::bit_vector> is(bv);
            c.decode(is, B.data(), n);
        }
        for (auto i = 0; i < n; i++) {
            ASSERT_EQ(B[i], A[i]);
        }
    }
}

TEST(bit_stream, fixed)
{
    size_t n = 20;
    std::mt19937 gen(4711);
    std::uniform_int_distribution<uint64_t> dis(1, 100);

    for (size_t i = 0; i < n; i++) {
        size_t len = dis(gen);
        std::vector<uint32_t> A(len);
        for (size_t j = 0; j < len; j++)
            A[j] = dis(gen);
        std::sort(A.begin(), A.end());
        auto last = std::unique(A.begin(), A.end());
        auto n = std::distance(A.begin(), last);
        coder::fixed<7> c;
        sdsl::bit_vector bv;
        {
            bit_ostream<sdsl::bit_vector> os(bv);
            c.encode(os, A.data(), n);
        }
        std::vector<uint32_t> B(n);
        {
            bit_istream<sdsl::bit_vector> is(bv);
            c.decode(is, B.data(), n);
        }
        for (auto i = 0; i < n; i++) {
            ASSERT_EQ(B[i], A[i]);
        }
    }
}

TEST(bit_stream, aligned_fixed)
{
    size_t n = 20;
    std::mt19937 gen(4711);
    std::uniform_int_distribution<uint64_t> dis(1, 100);

    for (size_t i = 0; i < n; i++) {
        size_t len = dis(gen);
        std::vector<uint32_t> A(len);
        for (size_t j = 0; j < len; j++)
            A[j] = dis(gen);
        std::sort(A.begin(), A.end());
        auto last = std::unique(A.begin(), A.end());
        auto n = std::distance(A.begin(), last);
        coder::aligned_fixed<uint32_t> c;
        sdsl::bit_vector bv;
        {
            bit_ostream<sdsl::bit_vector> os(bv);
            c.encode(os, A.data(), n);
        }
        std::vector<uint32_t> B(n);
        {
            bit_istream<sdsl::bit_vector> is(bv);
            c.decode(is, B.data(), n);
        }
        for (auto i = 0; i < n; i++) {
            ASSERT_EQ(B[i], A[i]);
        }
    }

    for (size_t i = 0; i < n; i++) {
        size_t len = dis(gen);
        std::vector<uint8_t> A(len);
        for (size_t j = 0; j < len; j++)
            A[j] = dis(gen);
        std::sort(A.begin(), A.end());
        auto last = std::unique(A.begin(), A.end());
        auto n = std::distance(A.begin(), last);
        coder::aligned_fixed<uint8_t> c;
        sdsl::bit_vector bv;
        {
            bit_ostream<sdsl::bit_vector> os(bv);
            c.encode(os, A.data(), n);
        }
        std::vector<uint8_t> B(n);
        {
            bit_istream<sdsl::bit_vector> is(bv);
            c.decode(is, B.data(), n);
        }
        for (auto i = 0; i < n; i++) {
            ASSERT_EQ(B[i], A[i]);
        }
    }
}

TEST(bit_stream, zlib)
{
    size_t n = 20;
    std::mt19937 gen(4711);
    std::uniform_int_distribution<uint64_t> dis(1, 100000);

    for (size_t i = 0; i < n; i++) {
        size_t len = dis(gen);
        std::vector<uint32_t> A(len);
        for (size_t j = 0; j < len; j++)
            A[j] = dis(gen);
        std::sort(A.begin(), A.end());
        auto last = std::unique(A.begin(), A.end());
        auto n = std::distance(A.begin(), last);
        coder::zlib<6> c;
        sdsl::bit_vector bv;
        {
            bit_ostream<sdsl::bit_vector> os(bv);
            c.encode(os, A.data(), n);
        }
        std::vector<uint32_t> B(n);
        {
            bit_istream<sdsl::bit_vector> is(bv);
            c.decode(is, B.data(), n);
        }
        for (auto i = 0; i < n; i++) {
            ASSERT_EQ(B[i], A[i]);
        }
    }
}

TEST(bit_stream, zlib_uint8)
{
    size_t n = 20;
    std::mt19937 gen(4711);
    std::uniform_int_distribution<uint64_t> dis(1, 200);

    for (size_t i = 0; i < n; i++) {
        size_t len = dis(gen);
        std::vector<uint8_t> A(len);
        for (size_t j = 0; j < len; j++)
            A[j] = dis(gen);
        std::sort(A.begin(), A.end());
        auto last = std::unique(A.begin(), A.end());
        auto n = std::distance(A.begin(), last);
        coder::zlib<6> c;
        sdsl::bit_vector bv;
        {
            bit_ostream<sdsl::bit_vector> os(bv);
            c.encode(os, A.data(), n);
        }
        std::vector<uint8_t> B(n);
        {
            bit_istream<sdsl::bit_vector> is(bv);
            c.decode(is, B.data(), n);
        }
        for (auto i = 0; i < n; i++) {
            ASSERT_EQ(B[i], A[i]);
        }
    }
}

// TEST(bit_stream, lzma)
// {
//     size_t n = 20;
//     std::mt19937 gen(4711);
//     std::uniform_int_distribution<uint64_t> dis(1, 100000);

//     for (size_t i = 0; i < n; i++) {
//         size_t len = dis(gen);
//         std::vector<uint32_t> A(len);
//         for (size_t j = 0; j < len; j++)
//             A[j] = dis(gen);
//         std::sort(A.begin(), A.end());
//         auto last = std::unique(A.begin(), A.end());
//         auto n = std::distance(A.begin(), last);
//         coder::lzma<2> c;
//         sdsl::bit_vector bv;
//         {
//             bit_ostream<sdsl::bit_vector> os(bv);
//             c.encode(os, A.data(), n);
//         }
//         std::vector<uint32_t> B(n);
//         {
//             bit_istream<sdsl::bit_vector> is(bv);
//             c.decode(is, B.data(), n);
//         }
//         for (auto i = 0; i < n; i++) {
//             ASSERT_EQ(B[i], A[i]);
//         }
//     }
// }

TEST(bit_stream, lz4_uint8)
{
    size_t n = 20;
    std::mt19937 gen(4711);
    std::uniform_int_distribution<uint64_t> dis(1, 200);

    for (size_t i = 0; i < n; i++) {
        size_t len = dis(gen);
        std::vector<uint8_t> A(len);
        for (size_t j = 0; j < len; j++)
            A[j] = dis(gen);
        std::sort(A.begin(), A.end());
        auto last = std::unique(A.begin(), A.end());
        auto n = std::distance(A.begin(), last);
        coder::lz4hc<9> c;
        sdsl::bit_vector bv;
        {
            bit_ostream<sdsl::bit_vector> os(bv);
            c.encode(os, A.data(), n);
        }
        std::vector<uint8_t> B(n);
        {
            bit_istream<sdsl::bit_vector> is(bv);
            c.decode(is, B.data(), n);
        }
        for (auto i = 0; i < n; i++) {
            ASSERT_EQ(B[i], A[i]);
        }
    }
}

TEST(bit_stream, bzip2)
{
    size_t n = 20;
    std::mt19937 gen(4711);
    std::uniform_int_distribution<uint64_t> dis(1, 200);

    for (size_t i = 0; i < n; i++) {
        size_t len = dis(gen);
        std::vector<uint8_t> A(len);
        for (size_t j = 0; j < len; j++)
            A[j] = dis(gen);
        std::sort(A.begin(), A.end());
        auto last = std::unique(A.begin(), A.end());
        auto n = std::distance(A.begin(), last);
        coder::bzip2<9> c;
        sdsl::bit_vector bv;
        {
            bit_ostream<sdsl::bit_vector> os(bv);
            c.encode(os, A.data(), n);
        }
        std::vector<uint8_t> B(n);
        {
            bit_istream<sdsl::bit_vector> is(bv);
            c.decode(is, B.data(), n);
        }
        for (auto i = 0; i < n; i++) {
            ASSERT_EQ(B[i], A[i]);
        }
    }
}

TEST(bit_stream, elias_fano)
{
    size_t n = 20;
    std::mt19937 gen(4711);
    std::uniform_int_distribution<uint64_t> dis(1, 2000);

    for (size_t i = 0; i < n; i++) {
        size_t len = dis(gen);
        std::vector<uint8_t> A(len);
        for (size_t j = 0; j < len; j++)
            A[j] = dis(gen);
        std::sort(A.begin(), A.end());
        auto last = std::unique(A.begin(), A.end());
        auto n = std::distance(A.begin(), last);
        coder::elias_fano c;
        sdsl::bit_vector bv;
        {
            bit_ostream<sdsl::bit_vector> os(bv);
            c.encode(os, A.data(), n, 2000);
        }
        std::vector<uint8_t> B(n);
        {
            bit_istream<sdsl::bit_vector> is(bv);
            c.decode(is, B.data(), n, 2000);
        }
        for (auto i = 0; i < n; i++) {
            ASSERT_EQ(B[i], A[i]);
        }
    }
}


TEST(factor_coder, factor_coder_blocked_subdict)
{
    const uint32_t factorization_blocksize = 64 * 1024;
    const uint64_t dict_size_bytes = 256 * 1024 * 1024;
    const uint64_t literal_threshold = 3;
    const uint64_t dict_segment_size_bytes = 1024;
    const uint64_t dict_pointer_width = utils::CLog2<dict_size_bytes>();
    const uint64_t dict_page_size = 16 * dict_segment_size_bytes;
    const uint64_t in_page_offset_width = utils::CLog2<dict_page_size>();
    const uint64_t page_ptr_width = dict_pointer_width - in_page_offset_width;
    const uint64_t num_pages_in_dict = dict_size_bytes / dict_page_size;
    block_factor_data bfd(factorization_blocksize);

    std::mt19937 gen(4711);
    std::uniform_int_distribution<uint64_t> num_fact_dist(100, 7000);
    std::lognormal_distribution<> fact_len_dist(3, 0.9);
    std::uniform_int_distribution<uint64_t> lit_sym_dist(1, 255);
    std::poisson_distribution<> offset_fact_dist(dict_size_bytes/2);

    using p0_coder = factor_coder_blocked_subdict<literal_threshold,dict_page_size,num_pages_in_dict,coder::fixed<8>, coder::elias_fano, coder::fixed<dict_pointer_width>, coder::vbyte >; 

    size_t n = 20;
    for(size_t i=0;i<n;i++) {
        // generate data
        bfd.reset();
        bfd.num_factors = num_fact_dist(gen);
        for(size_t j=0;j<bfd.num_factors;j++) {
            auto flen = std::round(fact_len_dist(gen));
            bfd.lengths[j] = flen;
            if(flen <= literal_threshold) {
                for(size_t k=0;k<flen;k++) {
                    auto sym = lit_sym_dist(gen);
                    bfd.literals[bfd.num_literals++] = sym;
                }
            } else {
                // offset case
                auto offset = offset_fact_dist(gen);
                bfd.offsets[bfd.num_offsets++] = offset;
            }
        }

        // encode it 
        p0_coder c;
        sdsl::bit_vector bv;
        {
            bit_ostream<sdsl::bit_vector> os(bv);
            c.encode_block(os,bfd);
        }
        block_factor_data bfd_recover(factorization_blocksize);
        {
            bit_istream<sdsl::bit_vector> is(bv);
            c.decode_block(is,bfd_recover,bfd.num_factors);
        }
    }
}

int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
