add_library(SF_fluxField STATIC
    SF_flux.cpp
)
target_include_directories(SF_fluxField PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(SF_fluxField PUBLIC SF_headers)
set_lib_output(SF_fluxField flux)
