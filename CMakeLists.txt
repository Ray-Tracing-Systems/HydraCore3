cmake_minimum_required(VERSION 3.24)


project(LiteScene CXX)

add_library(litescene_project_compile_options INTERFACE)
target_compile_features(litescene_project_compile_options INTERFACE cxx_std_17)

option(LITESCENE_ENABLE_GLTF "Enable support for .gltf format" OFF)
option(LITESCENE_ENABLE_LOGGING "Enable logging" ON)

set(LITESCENE_PROJECT_ROOT ${PROJECT_SOURCE_DIR})
set(LITESCENE_SRC_DIR ${LITESCENE_PROJECT_ROOT}/src)
set(LITESCENE_EXTERNAL_DIR ${LITESCENE_PROJECT_ROOT}/3rd_party)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})


#include(cmake/includes.cmake)

add_subdirectory(3rd_party)

add_subdirectory(${LITESCENE_SRC_DIR})
add_subdirectory(scn_mgr)