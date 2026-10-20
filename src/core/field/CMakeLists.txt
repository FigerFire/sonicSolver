add_library(SF_field STATIC
    SF_field.cpp
    SF_scalarField.cpp
)
target_include_directories(SF_field PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(SF_field PUBLIC
    ${CMAKE_SOURCE_DIR}/src/core/state
)
target_link_libraries(SF_field PUBLIC SF_headers SF_state SF_residual)
target_link_libraries(SF_field PUBLIC SF_fluxField)
set_lib_output(SF_field field)
