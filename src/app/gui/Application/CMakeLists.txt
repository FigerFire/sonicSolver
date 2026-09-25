set(SONIC_GUI_APP_ICON
    "${CMAKE_CURRENT_SOURCE_DIR}/../Resources/icons/sonicGui.icns"
)
set_source_files_properties(${SONIC_GUI_APP_ICON} PROPERTIES
    MACOSX_PACKAGE_LOCATION "Resources"
)

add_executable(sonicGui MACOSX_BUNDLE
    main.cpp
    SF_application.cpp
    SF_application.h
    ${SONIC_GUI_APP_ICON}
)

target_compile_definitions(sonicGui PRIVATE
    SONIC_GUI_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
)

target_link_libraries(sonicGui PRIVATE
    Qt${QT_VERSION_MAJOR}::Widgets
    SF_gui
    SF_guiSplash
)

if(TARGET sonicSolver)
    add_dependencies(sonicGui sonicSolver)
endif()

if(SONIC_GUI_VTK_FOUND AND COMMAND vtk_module_autoinit)
    vtk_module_autoinit(
        TARGETS sonicGui SF_guiVTKView
        MODULES ${VTK_LIBRARIES}
    )
endif()

set_target_properties(sonicGui PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    MACOSX_BUNDLE_BUNDLE_NAME "sonicFlow"
    MACOSX_BUNDLE_GUI_IDENTIFIER "com.sonicsolver.sonicFlow"
    MACOSX_BUNDLE_ICON_FILE "sonicGui.icns"
)
