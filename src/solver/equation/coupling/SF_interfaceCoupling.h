#pragma once

/// @file SF_interfaceCoupling.h
/// @brief 暴露 OneFluid interface transport、jump 与 source coupling contract。
///
/// Data flow:
///   level-set/interface state + flow field
///       -> interface transport/source/jump contribution
///       -> equation assembly and committed interface state
///
/// 调用者提供 runtime synchronization；本文件不拥有 timestep 或 runner authority。

#include "SF_interfaces.h"
#include "core/state/SF_state.h"

#include <memory>
#include <vector>

namespace SF {
class MultiBlockMesh;
namespace Physics::InterfaceModels { class Model; }
namespace Equation::Coupling {

/// @brief 将界面模型的 canonical 主标量登记到状态注册表。
void registerInterfaceState(
    Physics::InterfaceModels::Model& model,
    const Field& field,
    State::VariableRegistry& registry);

/// @brief 单 Field 的 OneFluid 界面方程生命周期适配器。
class InterfaceEquationCoupling final
    : public FDM::IEquationSystemCoupling {
public:
    InterfaceEquationCoupling(
        Physics::InterfaceModels::Model& model,
        State::VariableRegistry& variables,
        FDM::IExecutionRuntime& runtime);

    void beginStep(Field& field, double dt) override;
    void prepareRHS(Field& field, double dt) override;
    void assembleRHS(Field& field, Residual& residual, double dt) override;
    void preparePressureCorrection(Field& field) override;
    void commitStep(Field& field, double dt) override;
    State::VariableRegistry* variables(Field&) override;
    const FDM::ITransportModel* transportModel(
        const Field&) const override;
    const FDM::IInterfaceJumpCondition* interfaceJumpCondition(
        const Field&) const override;

private:
    void prepareInterfaceState(Field& field);

    Physics::InterfaceModels::Model& model_;
    State::VariableRegistry& variables_;
    FDM::IExecutionRuntime& runtime_;
};

/// @brief 多 patch OneFluid 界面方程生命周期和 canonical 标量同步。
class MultiPatchInterfaceEquationCoupling final
    : public FDM::IEquationSystemCoupling {
public:
    MultiPatchInterfaceEquationCoupling(
        MultiBlockMesh& mesh,
        const std::vector<int>& localPatchIds,
        FDM::IExecutionRuntime& runtime,
        std::vector<std::unique_ptr<Physics::InterfaceModels::Model>>& models,
        std::vector<State::VariableRegistry>& variables);

    void beginStep(const std::vector<Field*>& fields, double dt) override;
    void prepareRHS(const std::vector<Field*>& fields, double dt) override;
    void assembleRHS(const std::vector<Field*>& fields,
                     const std::vector<Residual*>& residuals, double dt) override;
    void preparePressureCorrection(
        const std::vector<Field*>& fields) override;
    void commitStep(const std::vector<Field*>& fields, double dt) override;
    State::VariableRegistry* variables(Field& field) override;
    const FDM::ITransportModel* transportModel(
        const Field& field) const override;
    const FDM::IInterfaceJumpCondition* interfaceJumpCondition(
        const Field& field) const override;

private:
    int phiHaloDepth() const;
    void preparePhiHalo(const char* operation);
    void publishCurvature(const char* operation);
    size_t findPatch(const Field& field) const;

    MultiBlockMesh& mesh_;
    const std::vector<int>& localPatchIds_;
    FDM::IExecutionRuntime& runtime_;
    std::vector<std::unique_ptr<Physics::InterfaceModels::Model>>& models_;
    std::vector<State::VariableRegistry>& variables_;
};

} // namespace Equation::Coupling
} // namespace SF
