add_library(SF_guiIO STATIC
    SF_configDocument.cpp
    SF_configDocument.h
    SF_project.cpp
    SF_project.h
    SF_IO.h
)

target_include_directories(SF_guiIO PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(SF_guiIO PUBLIC
    Qt${QT_VERSION_MAJOR}::Core
)

set_lib_output(SF_guiIO gui)
