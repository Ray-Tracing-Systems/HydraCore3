set(MODULE_NAME litescene)

set(MODULE_SOURCES
    cmesh4.cpp
    hydraxml.cpp
    scene.cpp
    loader.cpp
    geometry.cpp
)

set(MODULE_LIBS
    pugixml
)

set(MODULE_PUBLIC_LIBS
    LiteMath
)
