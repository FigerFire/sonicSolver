# Run 目标的 ownership：control-flow 解释器属于 run 层，而不是 system 层。
add_library(SF_run STATIC
    SF_planExecutor.cpp
    SF_planExecutor.h
)
target_include_directories(SF_run PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(SF_run PUBLIC
    SF_headers
)
set_lib_output(SF_run run)
