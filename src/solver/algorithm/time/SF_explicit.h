#pragma once

/// @file SF_explicit.h
/// @brief 显式时间离散的单 stage 数学与临时状态 contract。
///
/// Data flow:
///   State_n + dt + stageIndex -> RHS -> stage combination -> published state
///
/// 本文件不选择 governing equations、MPI topology 或 global timestep loop。

#include "SF_config.h"
#include "core/interfaces/SF_equationCoupling.h"
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

struct ScalarRKStorage {
    std::vector<double> q0, k1, k2, k3, k4;
};

struct RKStorage {
    std::vector<double> q0, k1, k2, k3, k4;
    std::vector<ScalarRKStorage> registered;
};

/// @brief Solver-owned storage spanning the stages of one explicit step.
struct Workspace {
    FDM::TimeRecipe recipe = FDM::builtInTimeRecipe(
        FDM::TimeRecipeId::ForwardEuler);
    bool active = false;
    int nextStage = 0;
    std::vector<RKStorage> patches;
};

/// @brief One forward-Euler update used by the explicit recipe and PISO predictor.
void forwardEuler(Field& field, const Residual& residual, double dt);

/// @brief Snapshot Q_n and registered variables required by the scheme.
void begin(Workspace& workspace,
           const std::vector<Field*>& fields,
           State::StateBundle& state,
           const FDM::TimeRecipe& recipe,
           FDM::IEquationSystemCoupling* equationSystem);

/// @brief Execute exactly one plan-selected explicit stage.
void executeStage(
    Workspace& workspace,
    int stageIndex,
    const std::vector<Field*>& fields,
    std::vector<PatchWorkspace>& workspaces,
    State::StateBundle& state,
    FDM::IEquationSystemCoupling* equationSystem,
    const AssembleRHS& assembleRHS,
    const Publish& publish,
    const Validate& validate);

} // namespace SF::Time::Explicit
