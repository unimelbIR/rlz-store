#pragma once

#include "utils.hpp"
#include "collection.hpp"

enum class remap_mode_t {
    only_reorder,
    reorder_preferred,
    reorder_preferred_hapax
};

template <uint32_t t_segment_size_bytes,uint32_t t_page_size_bytes,remap_mode_t t_mode>
class dict_segment_score {
public:
    static uint32_t segment_size() {
        return t_segment_size_bytes;
    }
    static uint32_t page_size() {
        return t_page_size_bytes;
    }
    static remap_mode_t remap_mode() {
        return t_mode;
    }
    static std::string remap_mode_str() {
        if(t_mode == remap_mode_t::only_reorder) return "-rom=or";
            
        if(t_mode == remap_mode_t::reorder_preferred) return "-rom=rp";
        return "-rom=rpx";
    }
    static std::string type()
    {
        return "dict_segment_score-s=" + std::to_string(t_segment_size_bytes) + "-p=" + std::to_string(t_page_size_bytes);
    }

    static std::string file_name(collection& col, uint64_t size_in_bytes)
    {
        auto size_in_mb = size_in_bytes / (1024 * 1024);
        return col.path + "/index/" + type() + "-" + std::to_string(size_in_mb) + ".sdsl";
    }

public:
    static void create(collection& , bool , size_t )
    {

    }

};
