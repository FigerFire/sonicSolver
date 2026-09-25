add_library(SF_solverInterfaces INTERFACE)
target_sources(SF_solverInterfaces INTERFACE
    SF_interfaces.h
    SF_observer.h
    SF_transportModel.h
    SF_equationCoupling.h
    SF_executionRuntime.h
    SF_boundaryPipeline.h
    SF_parallelCoordinatorContract.h
    SF_solverStepper.h
    SF_capability.h
    SF_immersedSystem.h
    SF_immersedConstraint.h
    SF_buildCapabilities.h
)

target_include_directories(SF_solverInterfaces INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(SF_solverInterfaces INTERFACE
    SF_headers
    SF_model
    SF_state
    SF_field
)
