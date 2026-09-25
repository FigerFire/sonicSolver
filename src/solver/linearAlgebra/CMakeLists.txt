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

# HYPRE 是 linear-algebra backend，因此它的实现、include 与 target ownership
# 都在本目录，而不是 infrastructure/mpi。infrastructure 只提供
# communicator/ownership；基础设施层不得反向依赖 solver。
find_package(MPI QUIET COMPONENTS CXX)
find_package(HYPRE QUIET CONFIG)
if(MPI_CXX_FOUND AND HYPRE_FOUND)
    target_sources(SF_linearAlgebra PRIVATE
        backend/hypre/SF_hypreBackend.cpp
    )
    target_include_directories(SF_linearAlgebra PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/backend/hypre
    )
    target_link_libraries(SF_linearAlgebra PRIVATE
        MPI::MPI_CXX
        HYPRE::HYPRE
    )
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
