#pragma once

/// @file SF_systemBuilder.h
/// @brief 从冻结配置构造 ResolvedSimulationSystem。

#include "SF_resolvedSimulationSystem.h"

namespace SF::System {

ResolvedSimulationSystem build(
    const FDM::SolverConfig& config,
    const BuildRequest& request);

} // namespace SF::System
