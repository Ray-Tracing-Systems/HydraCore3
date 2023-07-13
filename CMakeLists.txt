cmake_minimum_required (VERSION 3.8)
project (test)

option(USE_VULKAN "Enable Vulkan implementation" ON)

# get rid of annoying MSVC warnings.
add_definitions(-D_CRT_SECURE_NO_WARNINGS)

set (CMAKE_CXX_STANDARD 17)

include_directories(".")
include_directories("external")
if(WIN32)
  link_directories(${ADDITIONAL_LIBRARY_DIRS})
endif()

find_package(OpenMP)
if(USE_VULKAN)
  find_package(Vulkan)
  add_compile_definitions(USE_VOLK)
  add_compile_definitions(USE_VULKAN)
endif()

include_directories(${Vulkan_INCLUDE_DIR})
include_directories("../vkutils"
                    "../volk" 
                    "../LiteMath" 
                    "../loader_utils"
                    "../geo_utils"
                    "../RayTracing"
                    "../utils")
                    
link_directories("../volk")
link_directories("../RayTracing/lib")

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -Wno-attributes")
  message(STATUS "Debug is ON!")
  add_definitions(-D_DEBUG)
  #uncomment this to detect broken memory problems via gcc sanitizers
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fsanitize-address-use-after-scope -fno-omit-frame-pointer -fsanitize=leak -fsanitize=undefined -fsanitize=bounds-strict")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -fPIC -Wno-attributes")
endif()

set(CORE_SRC main.cpp
             integrator_rt.cpp
             integrator_pt1.cpp
             integrator_pt2.cpp
             integrator_pt_scene.cpp
             mi_materials.cpp 
             imageutils.cpp 
             external/miniz.c
             ../LiteMath/Image2d.cpp
             ../loader_utils/pugixml.cpp
             ../loader_utils/hydraxml.cpp
             ../RayTracing/EmbreeRT.cpp
             ../vkutils/geom/cmesh.cpp)

set(VULKAN_SRC integrator_pt1_generated.cpp
               integrator_pt1_generated_init.cpp
               integrator_pt1_generated_ds.cpp
               ../vkutils/vk_utils.cpp
               ../vkutils/vk_copy.cpp
               ../vkutils/vk_buffers.cpp
               ../vkutils/vk_images.cpp
               ../vkutils/vk_context.cpp
               ../vkutils/vk_alloc_simple.cpp
               ../vkutils/vk_pipeline.cpp
               ../vkutils/vk_descriptor_sets.cpp
               ../vkutils/ray_tracing/vk_rt_utils.cpp
               ../loader_utils/scene_mgr.cpp
               ../RayTracing/VulkanRTX.cpp
               ../vkutils/geom/vk_mesh.cpp)


if(USE_VULKAN)               
  add_executable(hydra ${CORE_SRC} ${VULKAN_SRC})
  set(LINKED_VULKAN ${Vulkan_LIBRARY})
else()
  add_executable(hydra ${CORE_SRC})
  set(LINKED_VULKAN "")
endif()


set_target_properties(hydra PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
target_link_libraries(hydra LINK_PUBLIC ${LINKED_VULKAN} volk dl   
                      embree3 embree_sse42 embree_avx embree_avx2 lexers simd sys tasking # embree shit
                      ) 
else()
target_link_libraries(hydra LINK_PUBLIC ${LINKED_VULKAN} OpenMP::OpenMP_CXX volk dl    
                      embree3 embree_sse42 embree_avx embree_avx2 lexers simd sys tasking # embree shit
                      ) 
endif()