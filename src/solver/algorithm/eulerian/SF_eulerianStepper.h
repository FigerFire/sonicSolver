#pragma once

/// @file SF_eulerianStepper.h
/// @brief Eulerian–Eulerian 的 SIMPLE/PISO/PIMPLE 压力耦合调度。

#include "SF_config.h"
#include "SF_phaseBoundary.h"
#include "solver/equation/eulerian/SF_equations.h"
#include "SF_phaseSource.h"
#include "core/interfaces/SF_solverStepper.h"
#include "solver/system/SF_runtimeRequirements.h"
#include "solver/equation/SF_assemblyPlan.h"
#include "solver/system/SF_stateRealizer.h"
#include "solver/run/SF_planExecutor.h"

#include <optional>

namespace SF::EulerianEulerian {

/// @brief 只调度耦合顺序；方程、边界、通信和临时量由独立模块负责。
class EulerianStepper final : public FDM::INavierStokesStepper {
public:
    /// 只依赖 WHAT / ORDER / RUNTIME；PISO 数学不变。
    /// @brief 只依赖 WHAT / HOW / ORDER / RUNTIME；PIMPLE 数学不变。
    /// @param numerics compiled HOW（time recipe + dt policy）。
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
    Equation::AssemblyPlanRegistry assemblyPlans_;
    PhaseEquationAssembler equations_;
    StepSummary lastSummary_;
    FDM::SolverServices services_;
    State::StateBundle* state_ = nullptr;
    bool servicesBound_ = false;
    bool solveStagesBound_ = false;
    bool pimpleStageActive_ = false;
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
