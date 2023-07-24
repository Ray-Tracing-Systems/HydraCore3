cmake_minimum_required (VERSION 3.15)
project (HydraCore3)

option(USE_VULKAN "Enable Vulkan implementation" ON)

# get rid of annoying MSVC warnings.
add_definitions(-D_CRT_SECURE_NO_WARNINGS)

set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP)

if(USE_VULKAN)

  if(CMAKE_SYSTEM_NAME STREQUAL Windows)
    set(VOLK_STATIC_DEFINES VK_USE_PLATFORM_WIN32_KHR)
  elseif(CMAKE_SYSTEM_NAME STREQUAL Linux)
    set(VOLK_STATIC_DEFINES VK_USE_PLATFORM_XLIB_KHR)
  elseif(CMAKE_SYSTEM_NAME STREQUAL Darwin)
    set(VOLK_STATIC_DEFINES VK_USE_PLATFORM_MACOS_MVK)
  endif()

  add_subdirectory(external/volk)

  add_compile_definitions(USE_VOLK)
  add_compile_definitions(USE_VULKAN)
endif()

include_directories(${Vulkan_INCLUDE_DIR})
include_directories("external"
                    "external/tinyexr"
                    "external/vkutils"
                    "external/volk"
                    "external/LiteMath"
                    "external/CrossRT")

link_directories("external/CrossRT/lib")

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -Wno-attributes")
  message(STATUS "Debug is ON!")
  add_definitions(-D_DEBUG)
  #uncomment this to detect broken memory problems via gcc sanitizers
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fsanitize-address-use-after-scope -fno-omit-frame-pointer -fsanitize=leak -fsanitize=undefined -fsanitize=bounds-strict")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -fPIC -Wno-attributes")
endif()


set(VULKAN_EXTERNAL_SRC
        external/vkutils/vk_utils.cpp
        external/vkutils/vk_copy.cpp
        external/vkutils/vk_buffers.cpp
        external/vkutils/vk_images.cpp
        external/vkutils/vk_context.cpp
        external/vkutils/vk_alloc_simple.cpp
        external/vkutils/vk_pipeline.cpp
        external/vkutils/vk_descriptor_sets.cpp
        external/vkutils/ray_tracing/vk_rt_utils.cpp
        external/vkutils/geom/vk_mesh.cpp
        external/CrossRT/VulkanRTX.cpp)

set(VULKAN_SRC
        integrator_pt1_generated.cpp
        integrator_pt1_generated_init.cpp
        integrator_pt1_generated_ds.cpp )

set(LOADER_EXTERNAL_SRC
        external/tinyexr/miniz.c
        external/LiteScene/scene_mgr.cpp
        external/LiteScene/pugixml.cpp
        external/LiteScene/hydraxml.cpp
        external/LiteScene/cmesh4.cpp
        external/vkutils/geom/cmesh.cpp)

set(OTHER_EXTERNAL_SRC
        external/LiteMath/Image2d.cpp
        external/CrossRT/EmbreeRT.cpp)

set(CORE_SRC
        main.cpp
        integrator_rt.cpp
        integrator_pt1.cpp
        integrator_pt2.cpp
        integrator_pt_scene.cpp
        mi_materials.cpp
        imageutils.cpp)


if(USE_VULKAN)               
  add_executable(hydra ${CORE_SRC}
                       ${VULKAN_EXTERNAL_SRC}
                       ${VULKAN_SRC}
                       ${LOADER_EXTERNAL_SRC}
                       ${OTHER_EXTERNAL_SRC})

  set(LINKED_VULKAN ${Vulkan_LIBRARY})
else()
  add_executable(hydra ${CORE_SRC}
                       ${LOADER_EXTERNAL_SRC}
                       ${OTHER_EXTERNAL_SRC})
  set(LINKED_VULKAN "")
endif()


set_target_properties(hydra PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")

set(EMBREE_LIBS embree3 embree_sse42 embree_avx embree_avx2 lexers simd sys tasking)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  target_link_libraries(hydra LINK_PUBLIC
                              ${LINKED_VULKAN} volk dl
                              ${EMBREE_LIBS}
                              )
else()
  target_link_libraries(hydra LINK_PUBLIC
                              ${LINKED_VULKAN} OpenMP::OpenMP_CXX volk dl
                              ${EMBREE_LIBS}
                              )
endif()