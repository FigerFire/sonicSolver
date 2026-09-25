#pragma once

/// @file SF_state.h
/// @brief 状态层统一公共接口。
///
/// 提供方程组变量布局、canonical 状态包以及 transported-variable
/// 注册表。模块外代码应包含本头文件，而非依赖具体 state 子组件。

#include "SF_stateLayout.h"
#include "SF_distributedField.h"
#include "SF_variableRegistry.h"
#include "SF_stateBundle.h"
