#pragma once

/// @file SF_presets.h
/// @brief COMPOSE — built-in 方程族 preset 入口（WHAT only）。
///
/// preset 只声明未知量/方程/约束；coupling、plan、backend、time recipe 由
/// 其它编译阶段决定。

#include "SF_buildRequest.h"
#include "SF_equationContribution.h"
#include "core/system/SF_equationIR.h"

#include <string>
#include <vector>

namespace SF::System::Compose {

/// @brief 单流体压力约束方程族（rho/U/rhoE + pressure-velocity constraint）。
void addPressureConstraintFluid(
    SystemCompositionBuilder& system, const PressureConstraintSpec& spec);

/// @brief rho = rho0 的连续方程 + 动量方程 + div(U)=0 约束。
void addConstantDensityFluid(
    SystemCompositionBuilder& system,
    const EquationCompositionConfig& composition, bool includeDiffusion);

/// @brief Eulerian-Eulerian 逐相方程 pack 与共享压力/体积分数约束。
void addEulerianEulerianTemplate(
    SystemCompositionBuilder& system, ResolvedSimulationSystem& resolved,
    const std::vector<std::string>& names);

} // namespace SF::System::Compose
