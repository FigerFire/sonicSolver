#pragma once

/// @file SF_turbulenceSystemContribution.h
/// @brief Turbulence-owned mathematical equation contribution contract.

#include <string>
#include <vector>

#include "core/system/SF_systemContribution.h"

namespace SF::Turbulence {

struct SystemContributionSpec {
    std::string model;
    std::vector<std::string> phases;
    bool eulerian = false;
    /// @brief 完整 PhaseSystem 声明顺序；selected phases 可以是子集或不同顺序。
    std::vector<std::string> phaseOrder;
};

void contribute(
    System::SystemContribution& system,
    const SystemContributionSpec& spec);

} // namespace SF::Turbulence
