enable_language(C)
enable_language(CXX)

if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
    set(CMAKE_POLICY_DEFAULT_CMP0135 NEW)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-O0 -g)
else()
    add_compile_options(-O2)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    add_compile_options(-g)
endif()

add_compile_options(-fPIC -w)

include(properties.cmake)
