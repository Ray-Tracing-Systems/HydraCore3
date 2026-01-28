enable_language(C)
enable_language(CXX)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-O0 -g)
else()
    add_compile_options(-O2)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    add_compile_options(-g)
endif()

if(NOT WIN32)
    add_compile_options(-fPIC)
endif()

include(properties.cmake)
