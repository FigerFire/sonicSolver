add_library(SF_executionRuntime STATIC
    SF_executionRuntime.cpp
    SF_localExecutionRuntime.cpp
)
target_include_directories(SF_executionRuntime PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src/core/interfaces
)
target_link_libraries(SF_executionRuntime PUBLIC
    SF_headers
    SF_solverInterfaces
    SF_state
)
set_lib_output(SF_executionRuntime execution)
