#pragma once

/// @file SF_explicit.h
/// @brief 通用显式时间积分 contract；给定 state、dt 与 RHS provider，执行
///        Euler、SSPRK3 或 RK4 的 stage snapshot、组合、发布与验证。
///
/// Data flow:
///   State_n + dt + RHS callback
///       -> stage state / RHS snapshot
///       -> explicit update and publish
///       -> State_{n+1}
///
/// 本文件不选择 governing equations、MPI topology 或 global timestep loop。

#include "SF_config.h"
#include "SF_interfaces.h"
#include "core/state/SF_state.h"
#include "solver/algorithm/SF_patchWorkspace.h"

#include <functional>
#include <vector>

namespace SF::Time::Explicit {

using SolverAlgorithm::PatchWorkspace;
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

} // namespace SF::Time::Explicit
