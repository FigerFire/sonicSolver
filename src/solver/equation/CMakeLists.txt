add_library(SF_equation STATIC
    SF_expression.cpp
    SF_assemblyPlan.cpp
    compressible/SF_compressible.cpp
    eulerian/SF_phase.cpp
    eulerian/SF_pressure.cpp
    eulerian/SF_timeStep.cpp
    eulerian/SF_transport.cpp
    eulerian/SF_turbulence.cpp
    eulerian/SF_workspace.cpp
)
target_include_directories(SF_equation PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/compressible
    ${CMAKE_CURRENT_SOURCE_DIR}/eulerian
)
target_link_libraries(SF_equation PUBLIC
    SF_headers
    SF_solverConfig
    SF_solverInterfaces
    SF_field
    SF_state
    SF_discretization
    SF_linearAlgebra
    SF_physics
    SF_turbulence
)
set_lib_output(SF_equation equation)
