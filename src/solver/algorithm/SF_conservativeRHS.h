#pragma once

/// @file SF_conservativeRHS.h
/// @brief Density-based single/multi patch 共用的 stage RHS 契约。

#include "SF_applicator.h"
#include "SF_config.h"
#include "core/interfaces/SF_solverStepper.h"
#include "core/state/SF_state.h"
#include "solver/equation/compressible/SF_compressible.h"
#include "solver/algorithm/SF_patchWorkspace.h"

#include <vector>

namespace SF::System { struct CompiledNumericalSystem; }

namespace SF::SolverAlgorithm::ConservativeRHS {

void validateTimestepState(
    const State::StateBundle& state,
    const FDM::SolverConfig& config,
    const System::CompiledNumericalSystem& numericalSystem);

void prepareBoundaryState(
    const std::vector<Field*>& fields,
    int conservativeHaloDepth,
    Boundary::Applicator& boundaryApplicator,
    FDM::SolverServices& services,
    double time,
    double dt);

void assembleAllPatches(
    const std::vector<Field*>& fields,
    std::vector<PatchWorkspace>& workspaces,
    const FDM::SolverConfig& config,
    const System::CompiledNumericalSystem& numericalSystem,
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

} // namespace SF::SolverAlgorithm::ConservativeRHS
