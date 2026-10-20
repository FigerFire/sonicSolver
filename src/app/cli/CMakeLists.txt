add_executable(sonicSolver
    main.cpp
    SF_cli.cpp
    SF_cli.h
    SF_resultCommands.cpp
    SF_resultCommands.h
)
target_link_libraries(sonicSolver PRIVATE SF_application)

set_target_properties(sonicSolver PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(sonicSolver PUBLIC OpenMP::OpenMP_CXX)
endif()
