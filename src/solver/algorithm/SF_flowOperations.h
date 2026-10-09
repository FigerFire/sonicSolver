#pragma once
#include "SF_field.h"
#include "solver/algorithm/time/SF_explicit.h"
#include "solver/algorithm/time/SF_stageGroup.h"
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
#include "solver/algorithm/pressure/SF_pressureOperators.h"
#include "solver/system/SF_runtimeRequirements.h"
#include "solver/system/SF_stateRealizer.h"

#include <functional>
#include <memory>
#include <string_view>
#include <utility>
#include <vector>

namespace SF {
namespace SolverAlgorithm {
/// Owns existing flow/RK/pressure numerical storage, never physical state or loops.
/// Registered operations retain this owner; services and compiled products are borrowed.
class FlowOperations : public std::enable_shared_from_this<FlowOperations> {
public:
    FlowOperations(FDM::SolverConfig,
        const System::ExecutableEquationSystem&, const System::CompiledNumericalSystem&,
        const System::CompiledSolvePlan&, const System::RuntimeRequirements&,
        State::StateBundle&, const FDM::SolverServices&, System::StateRealization&,
        const double& maximumTimeStep);
    void prepare();
    void bindOperations(Run::OpRegistry&);
    Time::TemporalCallbacks temporalParticipant();
    void bindStageSources();
    void prepareBoundaryState(const std::vector<Field*>& fields,double time,double dt);
    FDM::ImmersedConstraintResult& correctionResult() { return immersedCorrection_; }

private:
    FDM::SolverConfig config_;
    const System::ExecutableEquationSystem& executable_;
    const System::CompiledNumericalSystem& numerics_;
    std::unique_ptr<System::CompiledNumericalSystem> stageNumerics_;
    const System::CompiledNumericalSystem& rhsNumerics() const { return stageNumerics_?*stageNumerics_:numerics_; }
    const System::CompiledSolvePlan& solve_;
    const System::RuntimeRequirements& runtime_;
    Boundary::Applicator boundaryApplicator_;
    std::unique_ptr<SF::Equation::Compressible::System> equations_;
    std::unique_ptr<PressureBased::Corrector> conservativePressureCorrector_;
    std::unique_ptr<Pressure::PressureOperators> pressureOperators_;

    FDM::SolverServices services_;
    State::StateBundle* state_ = nullptr;
    std::vector<PatchWorkspace> workspaces_;
    bool solvePlanBound_ = false;
    std::vector<System::OpId> requiredOperations_;
    System::StateRealization& realizedState_;
    const double& maximumTimeStep_;
    /// @brief Stable handles into temporal backend arrays; no duplicated stage storage.
    std::vector<std::unique_ptr<State::DistributedFieldView>> temporalViews_;
    void bindTemporalStateViews();
    /// @brief solver-owned 显式 stage 存储（q0/k1..k4）；只按 layout 变化 resize。
    Time::Explicit::Workspace explicitWorkspace_;
    /// @brief 本步压力修正 summary；由 pressure 操作写入。
    PressureBased::CorrectionSummary pressureSummary_;
    FDM::ImmersedConstraintResult immersedCorrection_;

    void bindSolvePlan(const System::CompiledSolvePlan&);
    void resetStep();
    double proposeTimeStep() const;
    void beginExplicitStep();
    void executeExplicitStage(int);
    void bindState(State::StateBundle& state);
    void ensureWorkspaces(const std::vector<Field*>& fields);
    bool requiresOperation(std::string_view operation) const;
    void bindExplicitOps(Run::OpRegistry& operations,
                         const std::vector<Field*>& fields,
                         const double& maximumTimeStep,
                         Time::Explicit::Workspace& workspace);
    void bindPressureOps(Run::OpRegistry& operations,
                         Field& field,
                         const double& maximumTimeStep,
                         PressureBased::CorrectionSummary& summary);
    void bindConstantPressureOps(Run::OpRegistry& operations,
                                 const double& maximumTimeStep);

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
