#pragma once

/// @file SF_numerics.h
/// @brief numerics 模块总调度入口。
///
/// 本层只公开可供离散层选择的具体数值格式；方程项调度由顶层
/// `src/solver/discretization/SF_discretization.h` 负责。

#include "viscous/SF_viscous.h"
#include "methods/math/SF_newtonian.h"
#include "methods/numerics/state/SF_eulerState.h"
#include "immersed/SF_immersedOperators.h"
