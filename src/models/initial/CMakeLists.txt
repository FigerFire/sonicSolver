# Legacy 全局暂存配置（parser scratch state + makeSolverConfig 转换）已不再被
# 任何 target 编译：native case 解码现在直接写 typed CaseConfig。
# 文件保留在 app/application/model/legacy/ 供人工查证，见
# docs/runtime-authority-cleanup-report.md。

add_library(SF_initial STATIC SF_init.cpp)
target_include_directories(SF_initial PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(SF_initial PUBLIC
    SF_headers
    SF_state
    SF_solverConfig
    SF_field
)
set_lib_output(SF_initial initialCondition)
