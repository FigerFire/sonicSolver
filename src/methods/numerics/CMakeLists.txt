# ── 数值方法统一构建入口 ──

set(SF_NUMERICS_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/Flux
    ${CMAKE_CURRENT_SOURCE_DIR}/immersed
    ${CMAKE_CURRENT_SOURCE_DIR}/convection
    ${CMAKE_CURRENT_SOURCE_DIR}/structured
    ${CMAKE_CURRENT_SOURCE_DIR}/scalar
    ${CMAKE_CURRENT_SOURCE_DIR}/state
    ${CMAKE_CURRENT_SOURCE_DIR}/time
    ${CMAKE_CURRENT_SOURCE_DIR}/source
    ${CMAKE_CURRENT_SOURCE_DIR}/viscous
    ${CMAKE_SOURCE_DIR}/src/methods/math
    ${CMAKE_SOURCE_DIR}/src/methods/math/frame
    ${CMAKE_SOURCE_DIR}/src/methods/math/discrete
    ${CMAKE_SOURCE_DIR}/src/core/config
    ${CMAKE_SOURCE_DIR}/src/models/initial
    ${CMAKE_SOURCE_DIR}/src/models/ibm
    ${CMAKE_SOURCE_DIR}/src/models/physics
)

add_library(SF_riemann STATIC
    convection/SF_riemann.cpp
)
target_link_libraries(SF_riemann PUBLIC SF_headers SF_field)
set_lib_output(SF_riemann numerics)

add_library(SF_flux STATIC
    Flux/SF_roe.cpp
    Flux/SF_lw.cpp
    Flux/SF_lxF.cpp
    Flux/SF_sw.cpp
    Flux/SF_lowOrder.cpp
    Flux/SF_rusanovEOS.cpp
)
target_link_libraries(SF_flux PUBLIC SF_headers SF_field SF_riemann)
set_lib_output(SF_flux numerics)

add_library(SF_weno3 STATIC
    convection/SF_WENO3.cpp
)
target_link_libraries(SF_weno3 PUBLIC SF_headers SF_field SF_riemann SF_flux)
set_lib_output(SF_weno3 numerics)

add_library(SF_weno5 STATIC
    convection/SF_WENO5.cpp
)
target_link_libraries(SF_weno5 PUBLIC SF_headers SF_field SF_riemann SF_flux)
set_lib_output(SF_weno5 numerics)

add_library(SF_teno5 STATIC
    convection/SF_TENO5.cpp
)
target_link_libraries(SF_teno5 PUBLIC SF_headers SF_field SF_riemann SF_flux)
set_lib_output(SF_teno5 numerics)

add_library(SF_weno7 STATIC
    convection/SF_WENO7.cpp
)
target_link_libraries(SF_weno7 PUBLIC SF_headers SF_field SF_riemann SF_flux)
set_lib_output(SF_weno7 numerics)

foreach(target IN ITEMS
        SF_riemann SF_flux SF_weno3 SF_weno5 SF_teno5
        SF_weno7)
    target_include_directories(${target} PUBLIC ${SF_NUMERICS_INCLUDE_DIRS})
endforeach()
