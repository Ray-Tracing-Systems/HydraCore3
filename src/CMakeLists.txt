enable_language(CXX)
enable_language(C)

add_library(litescene_lib_compile_options INTERFACE)
target_link_libraries(litescene_lib_compile_options INTERFACE litescene_project_compile_options)
target_compile_options(litescene_lib_compile_options INTERFACE -m64)


if(LITESCENE_ENABLE_GLTF)
    target_compile_definitions(litescene_lib_compile_options INTERFACE LITESCENE_ENABLE_GLTF)
endif()
if(LITESCENE_ENABLE_LOGGING)
    target_compile_definitions(litescene_lib_compile_options INTERFACE LITESCENE_ENABLE_LOGGING)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(litescene_lib_compile_options INTERFACE -O0 -g)
else()
    target_compile_options(litescene_lib_compile_options INTERFACE -O1)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    target_compile_options(litescene_lib_compile_options INTERFACE -g)
endif()

add_library(hydra_xml INTERFACE)
target_sources(hydra_xml INTERFACE hydraxml.cpp)


include(properties.cmake)

add_library(${MODULE_NAME}
    ${MODULE_SOURCES}
)

target_include_directories(${MODULE_NAME} PUBLIC ${LITESCENE_PROJECT_ROOT}/include)

target_compile_options(${MODULE_NAME} PRIVATE -Wno-error=unused-variable)

target_link_libraries(${MODULE_NAME}
    PUBLIC ${MODULE_PUBLIC_LIBS} hydra_xml
    PRIVATE ${MODULE_LIBS} litescene_lib_compile_options 
)
