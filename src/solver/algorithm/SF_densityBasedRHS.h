#pragma once

/// @file SF_densityBasedRHS.h
/// @brief Density-based single/multi patch 共用的 stage RHS 契约。

#include "SF_applicator.h"
#include "SF_config.h"
#include "SF_interfaces.h"
#include "core/state/SF_state.h"
#include "solver/equation/compressible/SF_compressible.h"
#include "solver/algorithm/SF_patchWorkspace.h"

#include <vector>

namespace SF::SolverAlgorithm::DensityBasedRHS {

void validateTimestepState(
    const State::StateBundle& state,
    const FDM::SolverConfig& config);

void prepareBoundaryState(
    const std::vector<Field*>& fields,
    const FDM::SolverConfig& config,
    Boundary::Applicator& boundaryApplicator,
    FDM::SolverServices& services,
    double time,
    double dt);

void assembleAllPatches(
    const std::vector<Field*>& fields,
    std::vector<PatchWorkspace>& workspaces,
    const FDM::SolverConfig& config,
    Equation::Compressible::System& equations,
    State::StateBundle& state,
    FDM::SolverServices& services,
    Boundary::Applicator& boundaryApplicator,
    double stageTime);

void publishIntegratedState(
    const std::vector<Field*>& fields,
    const State::StateBundle& state,
    FDM::SolverServices& services);

void publishCorrectedConservativeState(
    bool wroteConservativeState,
    FDM::SolverServices& services);

void validateStateClosure(
    const std::vector<Field*>& fields,
    const State::StateBundle& state,
    FDM::SolverServices& services,
    const char* stage,
    bool labelPatches);

} // namespace SF::SolverAlgorithm::DensityBasedRHS
