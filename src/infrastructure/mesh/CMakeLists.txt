add_library(SF_mesh STATIC
    SF_meshConfig.h
    SF_LazyMeshLoader.cpp
    SF_mesh.cpp
    SF_cellFaceMesh.cpp
    communication/SF_communicationPlan.cpp
    MultiBlockMesh/SF_MultiBlockMesh.cpp
    decompose/SF_meshDecompose.cpp
    createMesh/SF_meshGen.cpp
    createMesh/SF_edges.cpp
    createMesh/SF_block.cpp
)
target_include_directories(SF_mesh PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/createMesh
    ${CMAKE_CURRENT_SOURCE_DIR}/decompose
    ${CMAKE_CURRENT_SOURCE_DIR}/MultiBlockMesh
    ${CMAKE_CURRENT_SOURCE_DIR}/communication
    ${CMAKE_CURRENT_SOURCE_DIR}/interface
    ${CMAKE_CURRENT_SOURCE_DIR}/slidingMesh
)
target_link_libraries(SF_mesh PUBLIC
    SF_headers
    SF_field
    SF_solverConfig
    SF_initial
    SF_io
    SF_physics
)
set_lib_output(SF_mesh mesh)
