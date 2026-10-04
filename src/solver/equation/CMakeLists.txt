add_library(SF_equation STATIC
    SF_expression.cpp
    compressible/SF_compressible.cpp
)
target_include_directories(SF_equation PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/compressible
)
target_link_libraries(SF_equation PUBLIC
    SF_headers
    SF_solverConfig
    SF_solverInterfaces
    SF_field
    SF_state
    SF_discretization
)
set_lib_output(SF_equation equation)

# Specialized Eulerian numerical kernels consume concrete phase closures;
# their supported mathematics are validated by the selected provider at compile time.
add_library(SF_eulerianEquationAdapter STATIC
    ../algorithm/eulerian/equation/SF_phase.cpp
    ../algorithm/eulerian/equation/SF_pressure.cpp
    ../algorithm/eulerian/equation/SF_timeStep.cpp
    ../algorithm/eulerian/equation/SF_transport.cpp
    ../algorithm/eulerian/equation/SF_workspace.cpp
)
target_link_libraries(SF_eulerianEquationAdapter PUBLIC
    SF_equation
    SF_physics
    SF_turbulence
    SF_mesh
    SF_linearAlgebra
)
set_lib_output(SF_eulerianEquationAdapter equation)
