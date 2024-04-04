#ifndef BIN_SERIALIZATION_H
#define BIN_SERIALIZATION_H
#include <ostream>
#include <istream>
#include "LiteMath.h"


namespace binary {

    bool is_little_endian();
    void serial_copy(const char *src, char *dst, unsigned size);
    void convert_to_native_order(const char *src, char *dst, unsigned size, bool from_big_endian);
    void convert_from_native_order(const char *src, char *dst, unsigned size, bool to_big_endian);

    template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>, void>>
    void write_ordered(std::ostream &dst, const T val, bool to_big_endian)
    {
        char buf[sizeof(T)];
        convert_from_native_order(reinterpret_cast<const char *>(&val), buf, sizeof(T), to_big_endian);
        dst.write(buf, sizeof(T));
    }

    template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>, void>>
    T read_ordered(std::istream &src, bool from_big_endian)
    {
        char buf[sizeof(T)];
        T val;
        src.read(buf, sizeof(T));
        convert_to_native_order(buf, reinterpret_cast<char *>(&val), sizeof(T), from_big_endian);
        return val;
    }

    template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>, void>>
    void write(std::ostream &dst, const T val)
    {
        char buf[sizeof(T)];
        serial_copy(reinterpret_cast<const char *>(&val), buf, sizeof(T));
        dst.write(buf, sizeof(T));
    }

    template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>, void>>
    T read(std::istream &src)
    {
        char buf[sizeof(T)];
        T val;
        src.read(buf, sizeof(T));
        serial_copy(buf, reinterpret_cast<char *>(&val), sizeof(T));
        return val;
    }

    template<typename T, int N>
    void write_vec(std::ostream &dst, const T &val)
    {
        for(int i = 0; i < N; ++i) {
            write<decltype(val.x)>(dst, val.M[i]);
        }
    }

    template<typename T, int N>
    T read_vec(std::istream &src)
    {
        T v{};
        for(int i = 0; i < N; ++i) {
            v.M[i] = read<decltype(v.x)>(src);
        }
        return v;
    }

}

#endif