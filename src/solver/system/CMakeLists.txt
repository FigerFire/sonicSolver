add_library(SF_solverSystem STATIC
    SF_systemBuilder.cpp
    SF_equationContribution.cpp
    SF_systemPrinter.cpp
    SF_systemValidator.cpp
    SF_stateRealizer.cpp
)

target_include_directories(SF_solverSystem PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(SF_solverSystem PUBLIC
    SF_headers
    SF_solverConfig
    SF_solverInterfaces
    SF_equation
)

set_lib_output(SF_solverSystem solverSystem)
