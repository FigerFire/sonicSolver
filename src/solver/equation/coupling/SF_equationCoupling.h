#pragma once

/// @file SF_equationCoupling.h
/// @brief 暴露 phase/model contribution 到 equation-system coupling 的 contract。
///
/// Data flow:
///   StateBundle / phase fields / model closures
///       -> transport/source/equation coupling
///       -> RHS and closure contribution
///
/// 调用者是 execution composition；本文件不拥有 timestep loop 或 runner selection。

#include "core/interfaces/SF_equationCoupling.h"
#include "core/interfaces/SF_executionRuntime.h"
#include "core/interfaces/SF_transportModel.h"
#include "core/state/SF_state.h"

#include <vector>

namespace SF {
class MultiBlockMesh;
namespace Physics::FluidStateModel { class HomogeneousPhaseChange; }
namespace Physics::Multiphase { class MultiPhaseModel; }
namespace Equation::Coupling {

void registerMixtureState(
    Physics::Multiphase::MultiPhaseModel& model,
    const Field& field,
    std::vector<double>& phaseMassRHS,
    State::VariableRegistry& registry);

class CompositeTransportProvider final : public FDM::ITransportModel {
public:
    void setTurbulence(FDM::ITransportModel* model);
    void setMultiPhase(Physics::Multiphase::MultiPhaseModel* model);
    void setInterfaceModel(FDM::ITransportModel* model);
    bool active() const;
    void applyBoundary(const Field& field) override;
    void correct(const Field& field, double dt) override;
    double dynamicViscosity(
        const Field& field, int i, int j, int k,
        double laminarMu) const override;
    std::vector<std::string> distributedReadFields() const override;
    std::vector<std::string> distributedWriteFields() const override;
    int distributedHaloDepth() const override;
private:
    FDM::ITransportModel* turbulence_ = nullptr;
    Physics::Multiphase::MultiPhaseModel* multiPhase_ = nullptr;
    FDM::ITransportModel* interfaceModel_ = nullptr;
};

class HomogeneousPhaseChangeProvider final
    : public FDM::IEquationSystemCoupling {
public:
    explicit HomogeneousPhaseChangeProvider(
        Physics::FluidStateModel::HomogeneousPhaseChange& model);
    void assembleRHS(Field& field, Residual& residual, double dt) override;
    void commitStep(Field&, double) override;
private:
    Physics::FluidStateModel::HomogeneousPhaseChange& model_;
};

class MixtureEquationProvider final
    : public FDM::IEquationSystemCoupling {
public:
    MixtureEquationProvider(
        Physics::Multiphase::MultiPhaseModel& model,
        std::vector<double>& rhs,
        State::VariableRegistry& variables,
        FDM::IExecutionRuntime& runtime);
    void beginStep(Field& field, double dt) override;
    void prepareRHS(Field& field, double dt) override;
    void assembleRHS(Field& field, Residual& residual, double dt) override;
    void commitStep(Field& field, double dt) override;
    State::VariableRegistry* variables(Field&) override;
private:
    Physics::Multiphase::MultiPhaseModel& model_;
    std::vector<double>& rhs_;
    State::VariableRegistry& variables_;
    FDM::IExecutionRuntime& runtime_;
};

class MultiPatchMixtureEquationProvider final
    : public FDM::IEquationSystemCoupling {
public:
    MultiPatchMixtureEquationProvider(
        MultiBlockMesh& mesh,
        const std::vector<int>& localPatchIds,
        FDM::IExecutionRuntime& runtime,
        std::vector<Physics::Multiphase::MultiPhaseModel>& models,
        std::vector<std::vector<double>>& rhs,
        std::vector<State::VariableRegistry>& variables,
        std::vector<CompositeTransportProvider>& transports);
    void beginStep(const std::vector<Field*>& fields, double dt) override;
    void prepareRHS(const std::vector<Field*>& fields, double dt) override;
    void assembleRHS(const std::vector<Field*>& fields,
                     const std::vector<Residual*>& residuals, double dt) override;
    void commitStep(const std::vector<Field*>& fields, double dt) override;
    State::VariableRegistry* variables(Field& field) override;
    const FDM::ITransportModel* transportModel(
        const Field& field) const override;
private:
    Execution::OperatorContract auxiliaryHaloContract(
        const char* operation) const;
    size_t findPatch(const Field& field) const;
    MultiBlockMesh& mesh_;
    const std::vector<int>& localPatchIds_;
    FDM::IExecutionRuntime& runtime_;
    std::vector<Physics::Multiphase::MultiPhaseModel>& models_;
    std::vector<std::vector<double>>& rhs_;
    std::vector<State::VariableRegistry>& variables_;
    std::vector<CompositeTransportProvider>& transports_;
};

} // namespace Equation::Coupling
} // namespace SF
