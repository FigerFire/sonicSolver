/// @file SF_singleFluidStepper.h
/// @brief 单流体可压缩算法编排器及其外部服务接口。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once
#include "SF_field.h"
#include "solver/algorithm/time/SF_explicit.h"
#include "solver/discretization/SF_discretization.h"
#include "solver/run/SF_planExecutor.h"
#include "SF_applicator.h"
#include "SF_config.h"
#include "core/interfaces/SF_solverStepper.h"
#include "SF_fluidStateModel.h"
#include "solver/equation/compressible/SF_compressible.h"
#include "core/state/SF_state.h"
#include "solver/algorithm/SF_patchWorkspace.h"
#include "solver/algorithm/pressureBased/SF_corrector.h"
#include "solver/algorithm/pressureBased/SF_kkt.h"
#include "solver/system/SF_runtimeRequirements.h"
#include "solver/system/SF_stateRealizer.h"

#include <functional>
#include <memory>
#include <string_view>
#include <utility>
#include <vector>

namespace SF {
namespace SolverAlgorithm {
/// @brief Navier-Stokes/Euler FDM step orchestrator.
///
/// `SingleFluidStepper` owns the common explicit predictor and optional pressure correction.
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
/// The compiled StageLoop invokes the Time::Explicit single-stage provider;
/// spatial assembly remains in equation/discretization providers.
///
/// The assembled conservative equation is:
/// `dQ/dt = -div(F_c) + div(F_v) + Source`.
class SingleFluidStepper : public FDM::INavierStokesStepper {
public:
    /// @brief Construct a solver from fully materialized FDM configuration.
    ///
/// This is the preferred constructor for new code and tests. The solver
/// copies the config and creates a default `Boundary::Applicator` from
/// `config.boundaries`.
    ///
    /// @param config Complete numerical/boundary/source configuration.
    /// 只依赖四类编译产物：WHAT / HOW / ORDER / RUNTIME。
    /// 不接收 ResolvedSimulationSystem，也不构造 AlgorithmContext。
    SingleFluidStepper(
        FDM::SolverConfig config,
        const System::ExecutableEquationSystem& equations,
        const System::CompiledNumericalSystem& numerics,
        const System::CompiledSolvePlan& solvePlan,
        const System::RuntimeRequirements& requirements);

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
    const System::ExecutableEquationSystem& executable_;
    const System::CompiledNumericalSystem& numerics_;
    const System::CompiledSolvePlan& solve_;
    const System::RuntimeRequirements& runtime_;
    Boundary::Applicator boundaryApplicator_;
    Equation::Compressible::System equations_;
    std::unique_ptr<PressureBased::Corrector> genericPisoCorrector_;
    std::unique_ptr<PressureBased::MonolithicKKT> monolithicKkt_;

    FDM::SolverServices services_;
    State::StateBundle* state_ = nullptr;
    std::vector<PatchWorkspace> workspaces_;
    bool servicesBound_ = false;
    bool solvePlanBound_ = false;
    std::vector<System::OpId> requiredOperations_;
    System::StateRealization realizedState_;
    /// @brief solver-owned operation registry；每步复用，不重新分配。
    Run::OpRegistry operations_;
    /// @brief solver-owned 显式 stage 存储（q0/k1..k4）；只按 layout 变化 resize。
    Time::Explicit::Workspace explicitWorkspace_;
    /// @brief 本步压力修正 summary；由 pressure 操作写入。
    PressureBased::CorrectionSummary pressureSummary_;
    bool correctionPerformed_ = false;
    bool correctionWroteConservative_ = false;
    std::string correctionDetail_;

    void bindState(State::StateBundle& state);
    void ensureWorkspaces(const std::vector<Field*>& fields);
    bool requiresOperation(std::string_view operation) const;
    void bindExplicitOps(Run::OpRegistry& operations,
                         const std::vector<Field*>& fields,
                         double maximumTimeStep,
                         Time::Explicit::Workspace& workspace);
    void bindPressureOps(Run::OpRegistry& operations,
                         Field& field,
                         double maximumTimeStep,
                         PressureBased::CorrectionSummary& summary);

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

} // namespace SolverAlgorithm
} // namespace SF
