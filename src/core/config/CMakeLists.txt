add_library(SF_solverConfig INTERFACE)
target_sources(SF_solverConfig INTERFACE
    SF_config.h
    SF_configTypes.h
    SF_equationComposition.h
    SF_numericsPolicy.h
    SF_runControl.h
)

target_include_directories(SF_solverConfig INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(SF_solverConfig INTERFACE
    SF_headers
    SF_state
)
