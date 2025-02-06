#ifndef SRC_LITESCENE_FORMAT_H_
#define SRC_LITESCENE_FORMAT_H_
#include <cstdio>
#include <string>
#include <iostream>

namespace ls::internal {
    
    template<typename ...Args>
    std::string format(const std::string &fmt, Args ...args)
    {
        unsigned long size = std::snprintf(nullptr, 0u, fmt.c_str(), args...) + 1;
        std::unique_ptr<char[]> ptr{new char[size]};
        std::snprintf(ptr.get(), size, fmt.c_str(), args...);
        return std::string(ptr.get());
    }
#ifdef LITESCENE_ENABLE_LOGGING
    template<typename ...Args>
    inline void log_error(const std::string &fmt, Args ...args)
    {
        std::cerr << "[LiteScene|Error] " << format(fmt, args...) << std::endl;
    }
#else
    template<typename ...Args>
    inline void log_error(const std::string &, Args ...) {}
#endif
}

#endif