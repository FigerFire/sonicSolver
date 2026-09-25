# GUI 热键静态库。
# 这里只保存键盘输入到业务模块的映射，不直接操作 MainWindow 内部控件。
add_library(SF_guiHotKey STATIC
    SF_hotKey.cpp
    SF_hotKey.h
)

target_include_directories(SF_guiHotKey PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(SF_guiHotKey PUBLIC
    Qt${QT_VERSION_MAJOR}::Widgets
    SF_guiIO
)

set_lib_output(SF_guiHotKey gui)
