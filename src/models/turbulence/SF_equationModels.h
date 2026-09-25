#pragma once

/// @file SF_equationModels.h
/// @brief Eulerian 湍流 EquationSystem 的模型闭式调度接口。

#include "SF_equationSystem.h"

namespace SF::Turbulence {

void prepareKEpsilon(
    const Physics::PhaseSystems::PhaseSystem& system,
    const FDM::TurbulenceConfig& config,
    PhaseEquationState& state);

void prepareKOmegaSST(
    const Physics::PhaseSystems::PhaseSystem& system,
    const FDM::TurbulenceConfig& config,
    PhaseEquationState& state);

void prepareSmagorinsky(
    const Physics::PhaseSystems::PhaseSystem& system,
    const FDM::TurbulenceConfig& config,
    PhaseEquationState& state);

} // namespace SF::Turbulence
