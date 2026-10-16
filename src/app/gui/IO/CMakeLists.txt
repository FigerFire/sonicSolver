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
    SF_application
    Qt${QT_VERSION_MAJOR}::Core
)

set_lib_output(SF_guiIO gui)

if(BUILD_TESTS)
    add_executable(test_nativeGuiProject ${CMAKE_SOURCE_DIR}/test/nativeArchitectureClosure/test_gui.cpp)
    target_link_libraries(test_nativeGuiProject PRIVATE SF_guiIO)
    target_compile_definitions(test_nativeGuiProject PRIVATE SF_TEST_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
endif()
