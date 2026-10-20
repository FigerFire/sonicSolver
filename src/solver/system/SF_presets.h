#pragma once

/// @file SF_presets.h
/// @brief COMPOSE — built-in 物理方程 bundle 贡献入口。
///
/// equation/state 与默认 occurrences/bindings 显式贡献；coupling topology、time recipe
/// 由独立请求决定。

#include "SF_buildRequest.h"
#include "SF_equationContribution.h"
#include "core/system/SF_equationIR.h"

#include <string>
#include <vector>

namespace SF::System::Compose {

/// @brief 单流体压力约束方程族（rho/U/rhoE + pressure-velocity constraint）。
void addPressureConstraintFluid(
    SystemCompositionBuilder& system, const PressureConstraintSpec& spec,
    const EquationCompositionConfig& composition);

/// @brief rho = rho0 的连续方程 + 动量方程 + div(U)=0 约束。
void addConstantDensityFluid(
    SystemCompositionBuilder& system,
    const EquationCompositionConfig& composition, bool includeDiffusion);

/// @brief Eulerian-Eulerian 逐相方程 pack 与共享压力/体积分数约束。
void addEulerianEulerianTemplate(
    SystemCompositionBuilder& system, ResolvedSimulationSystem& resolved,
    const std::vector<std::string>& names,const std::string& referencePhase);

} // namespace SF::System::Compose
