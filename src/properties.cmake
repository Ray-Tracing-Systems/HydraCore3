set(MODULE_NAME litescene)

set(MODULE_SOURCES
    cmesh4.cpp
    hydraxml.cpp
    scene.cpp
    loader.cpp
    loader_mat.cpp
    geometry.cpp
)

set(MODULE_LIBS
    pugixml
    tinyobjloader
)

set(MODULE_PUBLIC_LIBS
    LiteMath
)
