add_subdirectory(time)
add_subdirectory(workflow)

add_library(SF_solverAlgorithm STATIC
    SF_solverAlgorithm.h
    SF_compressible.cpp
    SF_compressible.h
    SF_densityBasedRHS.cpp
    SF_densityBasedRHS.h
    SF_densityBasedTime.cpp
    SF_densityBasedTime.h
    SF_solverAlgorithm.cpp
    densityBased/SF_correction.cpp
    densityBased/SF_densityBased.h
    pressureBased/SF_adapter.cpp
    pressureBased/SF_pressureBased.h
    pressureBased/SF_corrector.cpp
    pressureBased/SF_kkt.cpp
    pressureBased/SF_kktDofs.cpp
    pressureBased/SF_kktDofs.h
    pressureBased/SF_pressureJump.cpp
    pressureBased/eulerian/SF_pressureStepper.cpp
    immersed/SF_immersedStrategy.cpp
    immersed/SF_immersedStrategy.h
)
target_include_directories(SF_solverAlgorithm PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/densityBased
    ${CMAKE_CURRENT_SOURCE_DIR}/pressureBased
    ${CMAKE_CURRENT_SOURCE_DIR}/pressureBased/eulerian
    ${CMAKE_CURRENT_SOURCE_DIR}/time
    ${CMAKE_CURRENT_SOURCE_DIR}/workflow
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
)
set_lib_output(SF_solverAlgorithm solverAlgorithm)
