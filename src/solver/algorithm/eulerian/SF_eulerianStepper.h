#pragma once

/// @file SF_eulerianStepper.h
/// @brief Eulerian–Eulerian 的状态、数值 workspace 与 provider 回调。

#include "SF_config.h"
#include "SF_phaseBoundary.h"
#include "solver/algorithm/eulerian/equation/SF_equations.h"
#include "SF_phaseSource.h"
#include "core/interfaces/SF_solverStepper.h"
#include "solver/system/SF_runtimeRequirements.h"
#include "solver/system/SF_stateRealizer.h"
#include "solver/run/SF_planExecutor.h"

#include <optional>

namespace SF::EulerianEulerian {

/// @brief 绑定现有数值内核；不选择 equation 顺序或 coupling loop 次数。
class EulerianStepper final : public FDM::INavierStokesStepper {
public:
    /// @brief 只绑定已冻结的数学、数值选择与运行计划；PIMPLE 算术不变。
    /// @param numerics compiled WHICH（dt policy 与 phase numerical recipe）。
    EulerianStepper(Physics::PhaseSystems::PhaseSystem& phaseSystem,
                    FDM::SolverConfig config,
                    const System::ExecutableEquationSystem& equations,
                    const System::CompiledNumericalSystem& numerics,
                    const System::CompiledSolvePlan& solvePlan,
                    const System::RuntimeRequirements& requirements);

    /// @brief 把各相主状态、湍流量和压力 workspace 注册到 StateBundle。
    void registerState(State::StateBundle& state);

    void bindServices(FDM::SolverServices services) override;
    void bindSolvePlan(const System::CompiledSolvePlan& plan) override;
    void prepare(FDM::SolverState& state) override;
    FDM::StepResult advance(FDM::SolverState& state) override;
    /// @brief 返回由 compiled dt policy 决定的时间步上限。
    /// @param cfl compiled `CompiledNumericalSystem::dt.cfl`。
    double stableTimeStep(double cfl);
    const Turbulence::EquationSystem& turbulence() const {
        return turbulence_;
    }
    const StepSummary& lastSummary() const { return lastSummary_; }

private:
    Physics::PhaseSystems::PhaseSystem& system_;
    const System::ExecutableEquationSystem& executable_;
    const System::CompiledNumericalSystem& numerics_;
    const System::CompiledSolvePlan& solve_;
    const System::RuntimeRequirements& runtime_;
    FDM::SolverConfig config_;
    Turbulence::EquationSystem turbulence_;
    PhaseSolverWorkspace workspace_;
    Physics::PhaseSystems::PhaseBoundaryApplicator boundary_;
    Physics::PhaseSystems::PhaseSourceRegistry sourceRegistry_;
    PhaseEquationAssembler equations_;
    StepSummary lastSummary_;
    FDM::SolverServices services_;
    State::StateBundle* state_ = nullptr;
    bool servicesBound_ = false;
    bool solveStagesBound_ = false;
    /// @brief solver-owned operation registry；每步复用，不重新分配。
    Run::OpRegistry operations_;
    std::optional<LinearAlgebra::SolveResult> pendingPressureCorrection_;
    System::StateRealization realizedState_;

    void bindState(State::StateBundle& state);
    void registerOperations(
        Run::OpRegistry& operations,
        FDM::SolverState& solverState,
        double& dt);
    void applyBoundaryAndSynchronize();
    void synchronizeTurbulenceState();
    void synchronizePrimaryState();
    void synchronizePressureCorrection();
    void synchronizeMomentumDiagonal();
    void assembleCanonicalPhaseFlux();
};

} // namespace SF::EulerianEulerian
