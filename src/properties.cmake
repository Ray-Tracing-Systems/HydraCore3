set(MODULE_NAME LiteScene)

set(MODULE_SOURCES
    cmesh4.cpp
    #hydraxml.cpp
    scene.cpp
    loader.cpp
    loader_mat.cpp
    scene_2.cpp
    geometry.cpp
    mesh_load_obj.cpp
)

set(MODULE_LIBS

)

set(MODULE_PUBLIC_LIBS
    LiteMath
    tinyobjloader
    pugixml
)
