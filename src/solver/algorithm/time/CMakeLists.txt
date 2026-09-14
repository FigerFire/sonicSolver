add_library(SF_timeDriver STATIC
    SF_time.cpp
)

target_include_directories(SF_timeDriver PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(SF_timeDriver PUBLIC SF_headers SF_solverConfig)
set_lib_output(SF_timeDriver solverAlgorithm)
