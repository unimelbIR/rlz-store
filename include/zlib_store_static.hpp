#pragma once

#include "utils.hpp"
#include "collection.hpp"

#include "zlib.h"

using namespace std::chrono;

template <
uint32_t t_block_size = 2048,
uint32_t t_zlib_level = 6
>
class zlib_store_static {
private:
    sdsl::int_vector_mapper<8,std::ios_base::in> m_zlib_text;
    sdsl::int_vector<0>  m_offsets;
    sdsl::int_vector<0>  m_lens;
public:
    class builder;

    zlib_store_static() = delete;
    zlib_store_static(zlib_store_static&&) = default;
    zlib_store_static& operator=(zlib_store_static&&) = default;
    zlib_store_static(collection& col) 
    : m_zlib_text(col.file_map[KEY_ZLIB]) // (1) mmap factored text
    {
        // (2) load the block map
        LOG(INFO) << "\tLoad offsets";
        sdsl::load_from_file(m_offsets,col.file_map[KEY_ZLIB_OFFSETS]);
        // (3) load dictionary from disk
        LOG(INFO) << "\tLoad len";
        sdsl::load_from_file(m_lens,col.file_map[KEY_ZLIB_LEN]);
        LOG(INFO) << "len0 = " << m_lens[0];
        LOG(INFO) << "len1 = " << m_lens[1];
        LOG(INFO) << "len2 = " << m_lens[2];
        LOG(INFO) << "ZLib store ready";
    }

    std::vector<uint8_t>
    block(const size_t block_id) const {
        std::vector<uint8_t> block_content(t_block_size);
        z_stream strm;
	    strm.zalloc = Z_NULL;
	    strm.zfree = Z_NULL;
	    strm.opaque = Z_NULL;
	    strm.avail_in = 0;
	    strm.next_in = Z_NULL;
	    auto ret = inflateInit(&strm);
	    if (ret != Z_OK) {
	        LOG(FATAL) << "Unable to init zlib";
	    }

	    auto block_offset = m_offsets[block_id];
	    auto len = m_lens[block_id];
	    LOG(INFO) << "offset = " << block_offset;
	    LOG(INFO) << "len = " << len;
	    strm.avail_in = len;
	    strm.next_in = (uint8_t *)(m_zlib_text.data()+block_offset);
	    strm.avail_out = t_block_size;
	    strm.next_out = reinterpret_cast<uint8_t *>(block_content.data());
		inflate(&strm, Z_FINISH);
	    (void)inflateEnd(&strm);
        return block_content;
    }

    void
    block(const size_t block_id,std::vector<uint8_t>& block_content) const {

        return block_content;
    }
};

template <
uint32_t t_block_size = 2048,
uint32_t t_zlib_level = 6
>
class zlib_store_static<t_block_size,
                        t_zlib_level>::builder
{
    public:
        builder& set_rebuild(bool r) { rebuild = r; return *this; };
        builder& set_threads(uint8_t nt) { num_threads = nt; return *this; };

        template<class t_itr,class t_out>
        size_t zlib_compress_block(t_itr itr,t_out& out,size_t n) const
        {
            /* init zlib */
            z_stream strm;
		    strm.zalloc = Z_NULL;
		    strm.zfree = Z_NULL;
		    strm.opaque = Z_NULL;
		    auto ret = deflateInit(&strm, t_zlib_level);
		    if(ret != Z_OK) {
		    	LOG(FATAL) << "Unable to init zlib";
		    }
		    const uint32_t zlib_block_size = 4*t_block_size;
		    static uint8_t zlib_out_buf[zlib_block_size];
		    static uint8_t zlib_in_buf[t_block_size];
        	strm.avail_in = n;
        	std::copy(itr,itr+n,std::begin(zlib_in_buf));
        	strm.next_in = zlib_in_buf;
        	/* run until we are done processing */
        	size_t written_bytes = 0;
        	do {
        		strm.avail_out = zlib_block_size;
        		strm.next_out = zlib_out_buf;
        		ret = deflate(&strm, Z_FINISH);
        		auto have = zlib_block_size - strm.avail_out;
        		written_bytes += have;
        		/* output what we have */
        		for(size_t j=0;j<have;j++) {
        			out.push_back(zlib_out_buf[j]);
        		}
        	} while (strm.avail_out == 0);
            /* clean up */
            (void)deflateEnd(&strm);
        	return written_bytes;
        }

        zlib_store_static build_or_load(collection& col) const
        {
            auto start = hrclock::now();
            /* init input and output */
            {
                const sdsl::int_vector_mapper<8,std::ios_base::in> text(col.file_map[KEY_TEXT]);
                auto num_blocks = text.size() / t_block_size;
                auto zlib_file = col.path+"/index/"+KEY_ZLIB+"-"+std::to_string(t_block_size)+"-"+std::to_string(t_zlib_level)+".sdsl";
                auto zlib_out_stream = sdsl::write_out_buffer<8>::create(zlib_file);
                auto zlib_offset_file = col.path+"/index/"+KEY_ZLIB_OFFSETS+"-"+std::to_string(t_block_size)+"-"+std::to_string(t_zlib_level)+".sdsl";
                auto zlib_len_file = col.path+"/index/"+KEY_ZLIB_LEN+"-"+std::to_string(t_block_size)+"-"+std::to_string(t_zlib_level)+".sdsl";
                auto zlib_offset_stream = sdsl::write_out_buffer<0>::create(zlib_offset_file);
                auto zlib_len_out_stream = sdsl::write_out_buffer<0>::create(zlib_len_file);
                col.file_map[KEY_ZLIB] = zlib_file;
                col.file_map[KEY_ZLIB_OFFSETS] = zlib_offset_file;
                col.file_map[KEY_ZLIB_LEN] = zlib_len_file;

    		    /* process blocks */
    		    auto itr = text.begin();
    		    auto left = text.size() % t_block_size;
                for(size_t i=1;i<=num_blocks;i++) {
                	/* write offset */
                	zlib_offset_stream.push_back( zlib_out_stream.size() );
                	auto written_bytes = zlib_compress_block(itr,zlib_out_stream,t_block_size);
                	zlib_len_out_stream.push_back( written_bytes );
                	itr += t_block_size;
                }
                if(left) {
                	zlib_offset_stream.push_back( zlib_out_stream.size() );
                	auto written_bytes = zlib_compress_block(itr,zlib_out_stream,left);
                	zlib_len_out_stream.push_back( written_bytes );
                }
                sdsl::util::bit_compress(zlib_offset_stream);
                sdsl::util::bit_compress(zlib_len_out_stream);
            }

            auto stop = hrclock::now();
            LOG(INFO) << "ZLib construction complete. time = " << duration_cast<seconds>(stop - start).count() << " sec";

            return zlib_store_static(col);
        }
    private:
        bool rebuild = false;
        uint32_t num_threads = 1;
};