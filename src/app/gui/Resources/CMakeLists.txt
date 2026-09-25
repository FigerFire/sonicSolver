add_library(SF_guiResources STATIC
    SF_resources.qrc
)

target_link_libraries(SF_guiResources PUBLIC
    Qt${QT_VERSION_MAJOR}::Widgets
)

set_lib_output(SF_guiResources gui)
