/// @file SF_compressible.h
/// @brief 单流体可压缩算法编排器及其外部服务接口。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once
#include "SF_field.h"
#include "solver/discretization/SF_discretization.h"
#include "SF_applicator.h"
#include "SF_config.h"
#include "SF_interfaces.h"
#include "SF_equationSet.h"
#include "SF_flowAlgorithm.h"
#include "solver/equation/compressible/SF_compressible.h"
#include "core/state/SF_state.h"
#include "solver/algorithm/SF_patchWorkspace.h"
#include "solver/algorithm/pressureBased/SF_corrector.h"
#include "solver/system/SF_resolvedSimulationSystem.h"
#include "solver/system/SF_stateRealizer.h"

#include <functional>
#include <memory>
#include <utility>
#include <vector>

namespace SF::SolverAlgorithm {
/// @brief Navier-Stokes/Euler FDM step orchestrator.
///
/// `CompressibleAlgorithm` owns the common explicit predictor and optional pressure correction.
/// not own case parsing, mesh loading, MPI, IBM setup, or terminal output.
/// Those concerns are supplied as value configuration or small interfaces:
///
/// - `FDM::SolverConfig` selects schemes, transport constants, boundaries, and
///   source data.
/// - `Boundary::Applicator` 调用 `src/solver/boundary` 物理边界实现。
/// - `FDM::IExecutionRuntime` guarantees halo freshness at explicit stage boundaries.
/// - `FDM::IImmersedBoundary` enforces IBM ghost cells.
/// - `FDM::ISolverObserver` receives structured messages instead of the solver
///   printing directly.
///
/// Each numerical format remains a namespace operator (ddt/div/laplacian/Sp):
/// @code
///   RK4::ddt(field, dt, [&](Field& f) {
///       f.clearResidual();
///       WENO5::div(f, FDM::FluxSplitter::StegerWarming);
///       CENTRAL2::laplacian(f, mu, Pr);
///       SourceTerm::Sp(f, sourceConfig);
///   });
/// @endcode
///
/// The assembled conservative equation is:
/// `dQ/dt = -div(F_c) + div(F_v) + Source`.
class CompressibleAlgorithm : public FDM::INavierStokesStepper {
public:
    /// @brief Construct a solver from fully materialized FDM configuration.
    ///
/// This is the preferred constructor for new code and tests. The solver
/// copies the config and creates a default `Boundary::Applicator` from
/// `config.boundaries`.
    ///
    /// @param config Complete numerical/boundary/source configuration.
    CompressibleAlgorithm(
        FDM::SolverConfig config,
        const System::ResolvedSimulationSystem& system);

    void bindServices(FDM::SolverServices services) override;
    void bindSolvePlan(const System::CompiledSolvePlan& plan) override;
    void prepare(FDM::SolverState& state) override;

    /// @brief 推进 bundle 中全部 density-based patches。
    FDM::StepResult advance(FDM::SolverState& state) override;

    double physicalTime() const;

    double timeStep() const;

    /// @brief Access the immutable view of solver configuration.
    /// @return Current solver configuration value.
    const FDM::SolverConfig& config() const { return config_; }

private:
    FDM::SolverConfig config_;
    const System::ResolvedSimulationSystem& resolved_;
    Boundary::Applicator boundaryApplicator_;
    Equation::Compressible::System equations_;
    std::unique_ptr<FDM::IFlowAlgorithm> flowAlgorithm_;
    std::unique_ptr<PressureBased::Corrector> genericPisoCorrector_;

    FDM::SolverServices services_;
    State::StateBundle* state_ = nullptr;
    std::vector<PatchWorkspace> workspaces_;
    bool servicesBound_ = false;
    bool solveStagesBound_ = false;
    bool densitySolveBlockActive_ = false;
    bool pressurePredictorActive_ = false;
    bool pressureCorrectionActive_ = false;
    bool monolithicImmersedActive_ = false;
    System::OpId immersedPlanOperation_;
    System::StateRealization realizedState_;

    void bindState(State::StateBundle& state);
    void ensureWorkspaces(const std::vector<Field*>& fields);
    /// @brief 保持 pressure-based single-field route 的既有实现。
    void stepPressure(Field& field, double maximumTimeStep);
    void stepPressureGeneric(Field& field, double maximumTimeStep);
    void stepDensity(const std::vector<Field*>& fields,
                     double maximumTimeStep);

    /// @brief 通过注入管线或 legacy 适配顺序准备完整边界状态。
    void prepareBoundaryState(const std::vector<Field*>& fields,
                              double time, double dt);

    void correctTransportModel(const std::vector<Field*>& fields);
    /// @brief Emit a structured time-step message through the observer hook.
    void emitTimeStep();
    void emitConvectionContract();

    /// @brief Validate the conservative state after a time-integration stage.
    /// @param field Field to inspect.
    /// @param stage Human-readable stage label.
    void validateStateClosure(Field& field, const char* stage);
    void validateDensityStateClosure(const std::vector<Field*>& fields,
                                     const char* stage);
};

} // namespace SF::SolverAlgorithm
