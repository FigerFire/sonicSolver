#pragma once

/// @file SF_phaseMRF.h
/// @brief Eulerian 每相 MRF 方程贡献 provider。

#include "SF_valueTypes.h"
#include "SF_phaseSource.h"

#include <memory>
#include <vector>

namespace SF::Physics::PhaseSystems {

/// @brief 由显式旋转区域设置构造相方程 source provider。
std::unique_ptr<PhaseEquationSource> makePhaseMRFSource(
    std::vector<RotatingSetting> settings);

} // namespace SF::Physics::PhaseSystems
