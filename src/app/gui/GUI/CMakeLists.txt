# GUI 控件层静态库。
# 这里只放“显示控件”：主窗口、布尔开关、GEO/MESH 开关。
# Controller、IO、VTKView、Geometry 是其他子目录，各自有自己的 CMakeLists。
add_library(SF_guiWidgets STATIC
    SF_GUI.cpp
    SF_mainWindowLayout.cpp
    SF_mainWindowTree.cpp
    SF_GUI.h
    SF_mainWindow.h
    SF_switchButton.cpp
    SF_switchButton.h
    SF_viewModeSwitch.cpp
    SF_viewModeSwitch.h
)

target_include_directories(SF_guiWidgets PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# 链接关系：
# - Qt Widgets 提供 QLabel/QTreeWidget/QMenu 等控件。
# - SF_guiGeometry 提供中间几何编辑器。
# - SF_guiIO 提供 ConfigDocument/ProjectModule 等界面数据模型。
# - SF_guiResources 提供 qrc 图标资源。
# - SF_guiHotKey 提供 G/M/S 模块热键。
# - SF_guiVTKView 提供中间网格/结果显示窗口。
target_link_libraries(SF_guiWidgets PUBLIC
    Qt${QT_VERSION_MAJOR}::Widgets
    SF_guiGeometry
    SF_guiHotKey
    SF_guiIO
    SF_guiResources
    SF_guiVTKView
)

set_lib_output(SF_guiWidgets gui)
