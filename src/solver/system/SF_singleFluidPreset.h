#pragma once

/// @file SF_singleFluidPreset.h
/// @brief 内置守恒方程及默认四轴组合；不拥有数值执行。

#include "SF_buildRequest.h"
#include "SF_equationContribution.h"

namespace SF::System::Preset {

/// @brief 基础套餐的 typed 默认；用户 presence 与覆盖由控制层处理。
EquationCompositionConfig singleFluidDefaults();
/// @brief 直接构造守恒 NS 数学，返回原贡献容器；不分配 STATE 数组。
SystemContribution singleFluidEquations(const SingleFluidPresetSpec&);

/// @brief Install selected physical balances against explicitly selected solution targets.
///
/// WHAT selects balances; STATE supplies writable metadata and default targets.
/// Default HOW occurrences remain separate from the equation AST and its catalog.
void installSingleFluid(
    SystemCompositionBuilder& system,
    const SingleFluidPresetSpec& spec,
    const EquationCompositionConfig& composition);

} // namespace SF::System::Preset
