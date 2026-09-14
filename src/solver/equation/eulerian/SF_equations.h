#pragma once

/// @file SF_equations.h
/// @brief Eulerian 每相方程、Rhie–Chow 面通量和共享压力方程装配。

#include "solver/linearAlgebra/SF_linearAlgebra.h"
#include "solver/linearAlgebra/SF_distributedRowMap.h"
#include "solver/equation/eulerian/SF_workspace.h"
#include "solver/discretization/eulerian/SF_eulerian.h"
#include "SF_equationSystem.h"
#include "solver/equation/eulerian/SF_summary.h"

#include <memory>

namespace SF::EulerianEulerian {

/// @brief 数值方程装配层；物理闭式仍由 PhaseSystem 提供。
class PhaseEquationAssembler {
public:
    PhaseEquationAssembler(
        Physics::PhaseSystems::PhaseSystem& system,
        const FDM::SolverPropertiesConfig& config,
        PhaseSolverWorkspace& workspace);

    void setExecutionRuntime(FDM::IExecutionRuntime* runtime) {
        runtime_ = runtime;
    }
    void refreshRowMap();

    void initializeMomentumDiagonal(double dt);
    void buildMomentumInterpolatedFlux();
    void assembleContinuity(double dt, StepSummary& summary);
    void solveMomentumPredictors(double dt);
    void applySemiImplicitInterphase(double dt);
    LinearAlgebra::SolveResult solvePressureCorrection();
    void setPressureCorrection(const std::vector<double>& correction);
    void correctAllPhases();
    void correctCanonicalFaceFlux();
    void solvePhaseEnergy(double dt);
    /// @brief 求解已注册的湍流守恒输运方程块。
    void solveTurbulence(double dt);
    void attachTurbulence(Turbulence::EquationSystem* turbulence);
    void resetTurbulenceDiagnostics();
    int turbulenceIterations() const { return turbulenceIterations_; }
    double turbulenceResidual() const { return turbulenceResidual_; }
    LinearAlgebra::ReuseStatistics turbulenceReuseStatistics() const;

    /// @brief 由各相速度、声速和网格间距计算 Eulerian CFL 上限。
    double stableTimeStep(double cfl) const;
    /// @brief 返回相质量、相焓和相间强耦合源的时间尺度上限。
    double sourceTimeStep(double sourceCfl) const;
    /// @brief 实际推进前验证源项不会耗尽主状态。
    void validateSourceStep(double dt) const;
    const LinearAlgebra::ReuseStatistics& pressureReuseStatistics() const {
        return pressureSolver_.statistics();
    }

private:
    Physics::PhaseSystems::PhaseSystem& system_;
    const FDM::SolverPropertiesConfig& config_;
    PhaseSolverWorkspace& workspace_;
    LinearAlgebra::DistributedRowMap rowMap_;
    std::vector<std::unique_ptr<LinearAlgebra::SolverSession>>
        momentumSolvers_;
    std::vector<std::unique_ptr<LinearAlgebra::SolverSession>>
        energySolvers_;
    std::vector<std::unique_ptr<LinearAlgebra::SolverSession>>
        turbulenceSolvers_;
    Turbulence::EquationSystem* turbulence_ = nullptr;
    LinearAlgebra::SolverSession pressureSolver_;
    FDM::IExecutionRuntime* runtime_ = nullptr;
    double currentTimeStep_ = 0.0;
    int turbulenceIterations_ = 0;
    double turbulenceResidual_ = 0.0;

};

} // namespace SF::EulerianEulerian
