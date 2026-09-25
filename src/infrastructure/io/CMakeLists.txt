find_package(yaml-cpp REQUIRED)
add_library(SF_io STATIC
    serialization/SF_serialization.cpp
    registry/SF_IORegistry.cpp
    case/SF_case.cpp
    field/SF_field.cpp
    SF_meshDictReader.cpp
    private/SF_casePath.cpp
    private/SF_foamParser.cpp
    mesh/SF_sfmMesh.cpp
    mesh/SF_generatedMeshWriter.cpp
)
target_include_directories(SF_io PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)
target_include_directories(SF_io PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/mesh
    ${CMAKE_SOURCE_DIR}/src/infrastructure/mesh
    ${CMAKE_SOURCE_DIR}/src/infrastructure/mesh/createMesh
)
target_link_libraries(SF_io PUBLIC SF_headers SF_field SF_model yaml-cpp::yaml-cpp)
set_lib_output(SF_io io)

add_library(SF_resultWriter STATIC
    SF_resultWriter.h
    SF_vtkIndexWriter.cpp
    SF_vtkOutput.cpp
    SF_vtkPartitionWriter.cpp
    SF_vtkStructuredWriter.cpp
    private/SF_vtkUtils.h
)
target_include_directories(SF_resultWriter PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/src/infrastructure/mesh
)
target_link_libraries(SF_resultWriter PUBLIC
    SF_headers
    SF_field
    SF_solverInterfaces
)
set_lib_output(SF_resultWriter io)
