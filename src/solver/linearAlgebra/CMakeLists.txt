add_library(SF_linearAlgebra STATIC
    SF_globalDofSystem.cpp
    SF_linearAlgebra.cpp
    SF_distributedRowMap.cpp
    SF_distributedConstraintLayout.cpp
)
target_include_directories(SF_linearAlgebra PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/hypre
)
target_link_libraries(SF_linearAlgebra PUBLIC
    SF_headers
    SF_solverConfig
    SF_solverInterfaces
    SF_field
)

if(TARGET SF_hypreBackend)
    target_link_libraries(SF_linearAlgebra PRIVATE SF_hypreBackend)
    target_compile_definitions(SF_linearAlgebra PUBLIC SF_USE_HYPRE=1)
    message(STATUS "Pressure linear workflow enabled: HYPRE + MPI")
else()
    target_sources(SF_linearAlgebra PRIVATE hypre/SF_hypreUnavailable.cpp)
    target_compile_definitions(SF_linearAlgebra PUBLIC SF_USE_HYPRE=0)
    message(STATUS
        "Pressure linear workflow disabled: HYPRE/MPI not found; "
        "density workflows remain available")
endif()

set_lib_output(SF_linearAlgebra linearAlgebra)
