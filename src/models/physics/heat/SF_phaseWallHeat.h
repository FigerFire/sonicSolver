#pragma once

/// @file SF_phaseWallHeat.h
/// @brief Eulerian 壁面热通量方程贡献 provider。

#include "SF_valueTypes.h"
#include "SF_phaseSource.h"

#include <memory>
#include <vector>

namespace SF::Physics::PhaseSystems {

/// @brief 由壁面热通量设置构造相方程 source provider。
std::unique_ptr<PhaseEquationSource> makePhaseWallHeatSource(
    std::vector<WallHeatSetting> settings);

} // namespace SF::Physics::PhaseSystems
