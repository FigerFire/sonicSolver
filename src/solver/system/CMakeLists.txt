add_library(SF_solverSystem STATIC
    SF_systemBuilder.cpp
    SF_presets.cpp
    SF_presets.h
    SF_systemDescription.cpp
    SF_singleFluidPreset.cpp
    SF_stateRealization.cpp
    SF_stateRealization.h
    SF_pressureCoupling.cpp
    SF_pressureCoupling.h
    SF_couplingStatus.h
    SF_equationContribution.cpp
    SF_transformation.cpp
    SF_solvePlan.cpp
    SF_systemPrinter.cpp
    SF_systemValidator.cpp
    SF_stateRealizer.cpp
    SF_numericalCompiler.cpp
    SF_numericalCompiler.h
    SF_providerResolver.cpp
    SF_providerResolver.h
    SF_buildRequest.h
    SF_physicsTemplate.h
)

target_include_directories(SF_solverSystem PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# contribution 实现由各自 owning target 编译；system 只消费它们。
target_link_libraries(SF_solverSystem PUBLIC
    SF_headers
    SF_solverConfig
    SF_solverInterfaces
    SF_equation
    SF_physics
    SF_ibm
    SF_turbulence
)
target_link_libraries(SF_solverSystem PRIVATE
    SF_run
)

set_lib_output(SF_solverSystem solverSystem)
