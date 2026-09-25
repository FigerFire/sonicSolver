#pragma once

/// @file SF_phaseGravity.h
/// @brief Eulerian 每相重力方程贡献 provider。

#include "SF_valueTypes.h"
#include "SF_phaseSource.h"

#include <memory>
#include <vector>

namespace SF::Physics::PhaseSystems {

/// @brief 由显式重力设置构造相方程 source provider。
std::unique_ptr<PhaseEquationSource> makePhaseGravitySource(
    std::vector<ZoneVectorSetting> settings);

} // namespace SF::Physics::PhaseSystems
