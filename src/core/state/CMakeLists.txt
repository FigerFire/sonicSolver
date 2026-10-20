# 状态层为 header-only 契约模块；具体存储仍由 Field 所有。
add_library(SF_state INTERFACE)
target_sources(SF_state INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/SF_fluidStateModel.h
    ${CMAKE_CURRENT_SOURCE_DIR}/SF_idealGasDefaults.h
    ${CMAKE_CURRENT_SOURCE_DIR}/SF_state.h
    ${CMAKE_CURRENT_SOURCE_DIR}/SF_stateBundle.h
    ${CMAKE_CURRENT_SOURCE_DIR}/SF_stateLayout.h
    ${CMAKE_CURRENT_SOURCE_DIR}/SF_thermodynamicStateCache.h
    ${CMAKE_CURRENT_SOURCE_DIR}/SF_valueTypes.h
    ${CMAKE_CURRENT_SOURCE_DIR}/SF_variableRegistry.h
)
target_include_directories(SF_state INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(SF_state INTERFACE SF_headers)
