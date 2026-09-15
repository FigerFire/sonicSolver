add_library(SF_workflow STATIC
    SF_workflow.cpp
    SF_workflow.h
    SF_moduleGraph.cpp
    SF_moduleGraph.h
)

target_include_directories(SF_workflow PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(SF_workflow PUBLIC
    SF_headers
    SF_solverConfig
    SF_solverInterfaces
    SF_solverSystem
)

set_lib_output(SF_workflow algorithm/workflow)
