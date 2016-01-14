#pragma once

#include <sdsl/int_vector.hpp>
#include <string>

struct block_map_uncompressed {
    typedef typename sdsl::int_vector<>::size_type size_type;
    sdsl::int_vector<> m_block_offsets;
    sdsl::int_vector<> m_block_factors;
    sdsl::int_vector<> m_block_spep;

    static std::string type()
    {
        return "block_map_uncompressed";
    }

    block_map_uncompressed() = default;
    block_map_uncompressed(block_map_uncompressed&&) = default;
    block_map_uncompressed(block_map_uncompressed&) = default;
    block_map_uncompressed& operator=(const block_map_uncompressed&) = default;
    block_map_uncompressed& operator=(block_map_uncompressed&&) = default;

    block_map_uncompressed(collection& col)
    {
        LOG(INFO) << "\tLoad block offsets from file";
        sdsl::load_from_file(m_block_offsets, col.file_map[KEY_BLOCKOFFSETS]);
        if(utils::file_exists(col.file_map[KEY_BLOCKFACTORS])) 
            sdsl::load_from_file(m_block_factors, col.file_map[KEY_BLOCKFACTORS]);
        if(utils::file_exists(col.file_map[KEY_BLOCKSPEP])) 
            sdsl::load_from_file(m_block_spep, col.file_map[KEY_BLOCKSPEP]);
    }

    inline size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = NULL, std::string name = "") const
    {
        using namespace sdsl;
        structure_tree_node* child = structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += m_block_offsets.serialize(out, child, "offsets");
        written_bytes += m_block_factors.serialize(out, child, "num_factors");
        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    size_type size_in_bytes() const
    {
        return sdsl::size_in_bytes(*this);
    }

    inline size_type block_offset(size_t block_id) const
    {
        return m_block_offsets[block_id];
    }
    
    inline size_type block_factors(size_t block_id) const
    {
        return m_block_factors[block_id];
    }
    
    inline size_type block_spep(size_t block_id) const
    {
        return m_block_spep[block_id];
    }

    inline size_type num_blocks() const
    {
        return m_block_offsets.size();
    }
};
