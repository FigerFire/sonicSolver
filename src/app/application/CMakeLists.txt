add_library(SF_caseAdapter STATIC
    model/SF_model.cpp
    model/SF_output.cpp
    model/SF_nativeBinding.cpp
    model/compatibility/SF_compatibility.cpp
    model/compatibility/SF_documents.cpp
    model/compatibility/SF_nativeDecode.cpp
    model/compatibility/private/SF_sourceParserUtils.cpp
 )
target_include_directories(SF_caseAdapter PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/model/compatibility
    ${CMAKE_SOURCE_DIR}/src/infrastructure/io
    ${CMAKE_SOURCE_DIR}/src/infrastructure/mesh
    ${CMAKE_SOURCE_DIR}/src/infrastructure/mesh/createMesh
)
target_link_libraries(SF_caseAdapter PUBLIC SF_io SF_physics SF_solverConfig SF_timeDriver SF_resultWriter)
set_lib_output(SF_caseAdapter application)

add_library(SF_application STATIC
    SF_application.cpp
    SF_hostCapabilities.cpp
    SF_environment.cpp
    SF_inspection.cpp
    adapters/SF_ibmAdapters.cpp
    execution/SF_eulerian.cpp
    execution/SF_execution.cpp
    execution/SF_flowLoop.cpp
    execution/SF_multiPatch.cpp
    execution/SF_singleFluid.cpp
    model/SF_runtimeConfig.cpp
    output/SF_fields.cpp
    output/SF_report.cpp
)
target_include_directories(SF_application PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/adapters
    ${CMAKE_CURRENT_SOURCE_DIR}/execution
    ${CMAKE_CURRENT_SOURCE_DIR}/output
)
target_link_libraries(SF_application PUBLIC
    SF_headers
    SF_solverConfig
    SF_solverInterfaces
    SF_solverAlgorithm
    SF_solverSystem
    SF_timeDriver
    SF_state
    SF_field
    SF_boundary
    SF_io
    SF_caseAdapter
    SF_resultWriter
    SF_mpi
    SF_executionRuntime
    SF_mesh
    SF_ibm
    SF_physics
    SF_turbulence
)
set_lib_output(SF_application application)
