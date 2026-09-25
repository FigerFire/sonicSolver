add_library(SF_guiVTKView STATIC
    SF_widget.cpp
    SF_widget.h
    SF_VTKView.h
)

target_include_directories(SF_guiVTKView PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(SF_guiVTKView PUBLIC
    Qt${QT_VERSION_MAJOR}::Widgets
    SF_guiViewerControls
)

if(SONIC_GUI_VTK_FOUND)
    target_compile_definitions(SF_guiVTKView PUBLIC SONIC_GUI_WITH_VTK=1)
    target_link_libraries(SF_guiVTKView PUBLIC ${VTK_LIBRARIES})
endif()

set_lib_output(SF_guiVTKView gui)
