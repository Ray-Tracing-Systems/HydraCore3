
set(LITESCENE_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/3rd_party/pugixml.cpp
    ${CMAKE_CURRENT_LIST_DIR}/3rd_party/tinyexr/miniz.c
    ${CMAKE_CURRENT_LIST_DIR}/hydraxml.cpp
    ${CMAKE_CURRENT_LIST_DIR}/cmesh4.cpp
    ${CMAKE_CURRENT_LIST_DIR}/mesh_load_obj.cpp
    ${CMAKE_CURRENT_LIST_DIR}/scene.cpp
    ${CMAKE_CURRENT_LIST_DIR}/scene_mat.cpp
    ${CMAKE_CURRENT_LIST_DIR}/scene_tex.cpp
    ${CMAKE_CURRENT_LIST_DIR}/scene_convert.cpp
)

set(LITESCENE_VK_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/scene_mgr.cpp
)