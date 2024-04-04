#include "bin_serialization.h"
#include <iostream>
#include <algorithm>

namespace binary {

    bool is_little_endian()
    {
        int i = 1;
        return *(reinterpret_cast<char *>(&i));
    }

    void serial_copy(const char *src, char *dst, unsigned size)
    {
        std::copy(src, src + size, dst);
        if(is_little_endian()) {
            std::reverse(dst, dst + size);
        } 
    }

    void convert_to_native_order(const char *src, char *dst, unsigned size, bool from_big_endian)
    { 
        std::copy(src, src + size, dst);
        if(is_little_endian() ^ !from_big_endian) {
            std::reverse(dst, dst + size);
        } 
    }

    void convert_from_native_order(const char *src, char *dst, unsigned size, bool to_big_endian)
    {
        convert_to_native_order(src, dst, size, to_big_endian);
    }

}
