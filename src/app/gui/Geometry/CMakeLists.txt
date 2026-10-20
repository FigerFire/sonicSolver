add_library(SF_guiGeometry STATIC
    SF_editor.cpp
    SF_editor.h
    SF_Geometry.h
    SF_GeometryPanel.cpp
    SF_GeometryPanel.h
)

target_include_directories(SF_guiGeometry PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(SF_guiGeometry PUBLIC
    Qt${QT_VERSION_MAJOR}::Widgets
    SF_guiIO
)

set_lib_output(SF_guiGeometry gui)
