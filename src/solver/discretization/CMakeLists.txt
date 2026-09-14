add_library(SF_discretization INTERFACE)
target_sources(SF_discretization INTERFACE
    SF_discretization.h
    SF_deltaT.h
    time/SF_time.h
    convection/SF_convection.h
    convection/SF_WENO.h
    convection/SF_div.h
    diffusion/SF_diffusion.h
    diffusion/SF_laplacian.h
    source/SF_source.h
    source/SF_sourceTerm.h
    structured/SF_structured.h
    structured/SF_vectorCalculus.h
    eulerian/SF_eulerian.h
)
target_include_directories(SF_discretization INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/time
    ${CMAKE_CURRENT_SOURCE_DIR}/convection
    ${CMAKE_CURRENT_SOURCE_DIR}/diffusion
    ${CMAKE_CURRENT_SOURCE_DIR}/source
    ${CMAKE_CURRENT_SOURCE_DIR}/structured
    ${CMAKE_CURRENT_SOURCE_DIR}/eulerian
)
target_link_libraries(SF_discretization INTERFACE
    SF_headers
    SF_solverConfig
    SF_solverInterfaces
    SF_field
    SF_physics
    SF_riemann
    SF_flux
    SF_weno3
    SF_weno5
    SF_teno5
    SF_weno7
    SF_euler
    SF_rk4
)
