# 3D 预览控制静态库。
# 工具栏、显示模式枚举和 VTK 鼠标交互样式放在这里，VTKView 只负责装配。
add_library(SF_guiViewerControls STATIC
    SF_toolbar.cpp
    SF_toolbar.h
    SF_viewerControls.h
    SF_viewerTypes.h
)

target_include_directories(SF_guiViewerControls PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(SF_guiViewerControls PUBLIC
    Qt${QT_VERSION_MAJOR}::Widgets
)

if(SONIC_GUI_VTK_FOUND)
    target_compile_definitions(SF_guiViewerControls PUBLIC SONIC_GUI_WITH_VTK=1)
    target_sources(SF_guiViewerControls PRIVATE
        SF_interactorStyle.cpp
        SF_interactorStyle.h
    )
    target_link_libraries(SF_guiViewerControls PUBLIC ${VTK_LIBRARIES})
endif()

set_lib_output(SF_guiViewerControls gui)
