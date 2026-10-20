add_library(SF_guiController STATIC
    SF_controller.cpp
    SF_controller.h
)

target_include_directories(SF_guiController PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(SF_guiController PUBLIC
    Qt${QT_VERSION_MAJOR}::Core
    SF_guiIO
    SF_guiWidgets
)

set_lib_output(SF_guiController gui)
