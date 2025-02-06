enable_language(CXX)
enable_language(C)

add_library(scene_mgr INTERFACE)
target_sources(scene_mgr INTERFACE scene_mgr.cpp)

target_include_directories(scene_mgr INTERFACE .)