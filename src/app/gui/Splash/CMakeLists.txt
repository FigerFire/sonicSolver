add_library(SF_guiSplash STATIC
    SF_splash.cpp
    SF_splash.h
    SF_splash.qrc
)

target_include_directories(SF_guiSplash PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(SF_guiSplash PUBLIC
    Qt${QT_VERSION_MAJOR}::Widgets
    SF_guiResources
)

set_lib_output(SF_guiSplash gui)
