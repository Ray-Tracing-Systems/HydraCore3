include(properties.cmake)

add_library(${MODULE_NAME})

target_sources(${MODULE_NAME} PRIVATE 
	${MODULE_SOURCES}
)

target_include_directories(${MODULE_NAME} PUBLIC ./)

target_link_libraries(${MODULE_NAME}
	${MODULE_LIBS} litescene_project_compile_options
)
