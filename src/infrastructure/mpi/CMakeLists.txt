add_library(SF_mpi STATIC
    SF_globalDofResidual.cpp
    SF_haloExchange.cpp
    SF_parallelCoordinator.cpp
    SF_parallelContext.cpp
    backend/SF_mpiBackend.cpp
)
target_include_directories(SF_mpi PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/backend
    ${CMAKE_SOURCE_DIR}/src/infrastructure/mesh/communication
    ${CMAKE_SOURCE_DIR}/src/infrastructure/mesh/interface
    ${CMAKE_SOURCE_DIR}/src/core/interfaces
)
target_include_directories(SF_mpi PRIVATE
    ${CMAKE_SOURCE_DIR}/src/infrastructure/mesh
    ${CMAKE_SOURCE_DIR}/src/infrastructure/mesh/MultiBlockMesh
    ${CMAKE_SOURCE_DIR}/src/infrastructure/io
    ${CMAKE_SOURCE_DIR}/src/models/physics/multiphase
    ${CMAKE_SOURCE_DIR}/src/models/physics/phaseChange
    ${CMAKE_SOURCE_DIR}/src/models/ibm
    ${CMAKE_SOURCE_DIR}/src/core/config
    ${CMAKE_SOURCE_DIR}/src/methods/numerics
    ${CMAKE_SOURCE_DIR}/src/methods/numerics/Flux
    ${CMAKE_SOURCE_DIR}/src/methods/numerics/structured
)

target_link_libraries(SF_mpi PUBLIC SF_headers SF_field)

find_package(MPI QUIET COMPONENTS CXX)
if(MPI_CXX_FOUND)
    target_link_libraries(SF_mpi PRIVATE MPI::MPI_CXX)
    target_compile_definitions(SF_mpi PRIVATE SF_USE_MPI=1)
    message(STATUS "MPI enabled: ${MPI_CXX_COMPILER}")
else()
    target_compile_definitions(SF_mpi PRIVATE SF_USE_MPI=0)
    message(WARNING "MPI not found — SF_mpi builds as a serial stub.")
endif()

set_lib_output(SF_mpi mpi)
