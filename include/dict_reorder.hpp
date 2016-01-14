#pragma once

#include "utils.hpp"
#include "collection.hpp"

template <uint32_t t_page_size_bytes>
class dict_reorder {
public:
    static uint32_t page_size() {
        return t_page_size_bytes;
    }
    static std::string type()
    {
        return "dict_reorder-" + std::to_string(t_page_size_bytes);
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
