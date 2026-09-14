add_library(SF_boundary STATIC
    SF_applicator.h
    SF_pipeline.cpp
    SF_pipeline.h
    SF_boundary.cpp
    SF_thermodynamicClosure.cpp
    reconstruction/SF_reconstruction.cpp
    reconstruction/linear/SF_linear.cpp
    reconstruction/linear/SF_fixedValue.cpp
    reconstruction/linear/SF_zeroGradient.cpp
    reconstruction/linear/SF_symmetry.cpp
    reconstruction/linear/SF_empty.cpp
    reconstruction/ILW/SF_ILW.cpp
    reconstruction/ILW/SF_boundaryClosure.cpp
    reconstruction/ILW/SF_fixedValue.cpp
    reconstruction/ILW/SF_zeroGradient.cpp
    reconstruction/ILW/SF_symmetry.cpp
    reconstruction/ILW/SF_empty.cpp
)
target_include_directories(SF_boundary PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/law
        ${CMAKE_CURRENT_SOURCE_DIR}/reconstruction
        ${CMAKE_CURRENT_SOURCE_DIR}/reconstruction/linear
        ${CMAKE_CURRENT_SOURCE_DIR}/reconstruction/ILW
    ${CMAKE_SOURCE_DIR}/src/methods/math
    ${CMAKE_SOURCE_DIR}/src/methods/math/frame
    ${CMAKE_SOURCE_DIR}/src/methods/math/discrete
    ${CMAKE_SOURCE_DIR}/src/methods/numerics/state
    ${CMAKE_SOURCE_DIR}/src/solver/discretization
    ${CMAKE_SOURCE_DIR}/src/solver/discretization/structured
    ${CMAKE_SOURCE_DIR}/src/methods/numerics/structured
    ${CMAKE_SOURCE_DIR}/src/methods/numerics
    ${CMAKE_SOURCE_DIR}/src/models/ibm
)

target_link_libraries(SF_boundary PUBLIC
    SF_headers
    SF_field
    SF_solverConfig
    SF_solverInterfaces
)
set_lib_output(SF_boundary boundary)
