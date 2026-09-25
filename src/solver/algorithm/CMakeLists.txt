add_subdirectory(time)

add_library(SF_solverAlgorithm STATIC
    SF_singleFluidStepper.cpp
    SF_singleFluidStepper.h
    SF_conservativeRHS.cpp
    SF_conservativeRHS.h
    time/SF_explicit.cpp
    time/SF_explicit.h
    pressureBased/SF_corrector.cpp
    pressureBased/SF_kkt.cpp
    pressureBased/SF_kktDofs.cpp
    pressureBased/SF_kktDofs.h
    pressureBased/SF_pressureJump.cpp
    eulerian/SF_eulerianStepper.cpp
    eulerian/SF_eulerianStepper.h
    immersed/SF_immersedStrategy.cpp
    immersed/SF_immersedStrategy.h
    immersed/SF_constraintOps.cpp
    immersed/SF_constraintOps.h
)
target_include_directories(SF_solverAlgorithm PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/pressureBased
    ${CMAKE_CURRENT_SOURCE_DIR}/eulerian
    ${CMAKE_CURRENT_SOURCE_DIR}/time
    ${CMAKE_CURRENT_SOURCE_DIR}/immersed
    ${CMAKE_SOURCE_DIR}/src/models/ibm
)
target_link_libraries(SF_solverAlgorithm PUBLIC
    SF_headers
    SF_solverConfig
    SF_solverInterfaces
    SF_state
    SF_field
    SF_boundary
    SF_discretization
    SF_equation
    SF_linearAlgebra
    SF_physics
    SF_turbulence
    SF_run
)
set_lib_output(SF_solverAlgorithm solverAlgorithm)
