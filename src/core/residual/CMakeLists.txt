add_library(SF_residual STATIC
    SF_residual.cpp
)
target_include_directories(SF_residual PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(SF_residual PUBLIC SF_headers)
set_lib_output(SF_residual residual)
