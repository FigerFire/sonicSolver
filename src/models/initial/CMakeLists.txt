# Native case 解码直接写 typed CaseConfig；初始场构造不持有 parser 暂存状态。

add_library(SF_initial STATIC SF_init.cpp)
target_include_directories(SF_initial PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(SF_initial PUBLIC
    SF_headers
    SF_state
    SF_solverConfig
    SF_field
)
set_lib_output(SF_initial initialCondition)
