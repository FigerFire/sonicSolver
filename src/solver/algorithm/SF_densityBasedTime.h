#pragma once

/// @file SF_densityBasedTime.h
/// @brief Density-based 显式 stage algebra；空间 RHS 由调用者提供。

#include "SF_config.h"
#include "SF_interfaces.h"
#include "core/state/SF_state.h"
#include "solver/algorithm/SF_patchWorkspace.h"

#include <functional>
#include <vector>

namespace SF::SolverAlgorithm::DensityBasedTime {

using AssembleRHS = std::function<void(const std::vector<Field*>&,
                                       std::vector<PatchWorkspace>&, double)>;
using Publish = std::function<void(const std::vector<Field*>&)>;
using Validate = std::function<void(const std::vector<Field*>&, const char*)>;

/// Qn belongs to StateBundle.  q0/k workspaces and Field stage states are
/// temporary; StateBundle::time is never advanced inside this function.
void advance(
    const std::vector<Field*>& fields,
    std::vector<PatchWorkspace>& workspaces,
    State::StateBundle& state,
    FDM::TimeScheme scheme,
    FDM::IEquationSystemCoupling* equationSystem,
    const AssembleRHS& assembleRHS,
    const Publish& publish,
    const Validate& validate);

} // namespace SF::SolverAlgorithm::DensityBasedTime
