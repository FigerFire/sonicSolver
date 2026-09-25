add_library(SF_turbulence STATIC
    SF_turbulence.cpp
    SF_turbulenceSystemContribution.cpp
    SF_turbulenceSystemContribution.h
    SF_equationSystem.cpp
    SF_equationModelOps.cpp
    SF_kEpsilonEquation.cpp
    SF_kOmegaSSTEquation.cpp
    SF_smagorinskyEquation.cpp
    RAS/SF_RAS.cpp
    RAS/SF_kEpsilon.cpp
    RAS/SF_kOmegaSST.cpp
    LES/SF_LES.cpp
    DNS/SF_DNS.cpp
)
target_include_directories(SF_turbulence PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/RAS
    ${CMAKE_CURRENT_SOURCE_DIR}/LES
    ${CMAKE_CURRENT_SOURCE_DIR}/DNS
    ${CMAKE_SOURCE_DIR}/src/methods/numerics/viscous
    ${CMAKE_SOURCE_DIR}/src/solver/linearAlgebra
    ${CMAKE_SOURCE_DIR}/src/solver/discretization/structured
)
target_link_libraries(SF_turbulence PUBLIC
    SF_equation
    SF_headers
    SF_field
    SF_solverConfig
    SF_boundary
    SF_physics
    SF_linearAlgebra
)
set_lib_output(SF_turbulence turbulence)
