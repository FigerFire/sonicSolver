#pragma once

/// @file SF_pressureStepper.h
/// @brief Eulerian–Eulerian 的 SIMPLE/PISO/PIMPLE 压力耦合调度。

#include "SF_config.h"
#include "SF_phaseBoundary.h"
#include "solver/equation/eulerian/SF_equations.h"
#include "SF_phaseSource.h"
#include "SF_interfaces.h"
#include "solver/system/SF_resolvedSimulationSystem.h"
#include "solver/equation/SF_assemblyPlan.h"
#include "solver/system/SF_stateRealizer.h"

namespace SF::EulerianEulerian {

/// @brief 只调度耦合顺序；方程、边界、通信和临时量由独立模块负责。
class PressureStepper final : public FDM::INavierStokesStepper {
public:
    PressureStepper(Physics::PhaseSystems::PhaseSystem& phaseSystem,
                    FDM::SolverConfig config,
                    const System::ResolvedSimulationSystem& resolved);

    /// @brief 把各相主状态、湍流量和压力 workspace 注册到 StateBundle。
    void registerState(State::StateBundle& state);

    void bindServices(FDM::SolverServices services) override;
    void bindSolveStages(
        const std::vector<FDM::SolveStage>& stages) override;
    void prepare(FDM::SolverState& state) override;
    FDM::StepResult advance(FDM::SolverState& state) override;
    /// @brief 返回由 `[numerics].CFL/maxCo` 决定的时间步上限。
    double stableTimeStep(double cfl);
    const Turbulence::EquationSystem& turbulence() const {
        return turbulence_;
    }
    const StepSummary& lastSummary() const { return lastSummary_; }

private:
    Physics::PhaseSystems::PhaseSystem& system_;
    const System::ResolvedSimulationSystem& resolved_;
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
    bool turbulenceStageActive_ = false;
    System::StateRealization realizedState_;

    void bindState(State::StateBundle& state);
    StepSummary stepImpl(double dt);
    void applyBoundaryAndSynchronize();
    void synchronizeTurbulenceState();
    void synchronizePrimaryState();
    void synchronizePressureCorrection();
    void synchronizeMomentumDiagonal();
    void assembleCanonicalPhaseFlux();
};

} // namespace SF::EulerianEulerian
