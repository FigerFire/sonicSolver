#pragma once

/// @file SF_equationContribution.h
/// @brief Builtin、model 与 user 共用的 equation-system composition contract。

#include "SF_resolvedSimulationSystem.h"
#include "core/system/SF_systemContribution.h"

#include <memory>
#include <vector>

namespace SF::System {

/// @brief 向 raw system 注册普通数学对象，不执行 transformation 或 timestep。
class SystemCompositionBuilder {
public:
    SystemCompositionBuilder(
        RawEquationSystem& system,
        std::vector<TransformationDescriptor>& transformations,
        std::vector<LegacyExecutionPolicy>& policies,
        Provenance origin);

    void recordContribution(std::string id, std::string name);
    void addState(StateSymbol unknown);
    /// A request activates catalog metadata at freeze; it adds no equations/storage.
    void requireState(std::string id) { requiredStates.push_back(std::move(id)); }
    std::vector<std::string> requiredStates;
    void addEquation(
        EquationDescriptor descriptor, SF::Equation::Definition definition);
    void addEquation(
        EquationDescriptor descriptor, SF::Equation::Definition definition,
        Equation formula);
    void addEquation(Equation equation);
    void extendEquation(const std::string& equationId, SF::Equation::Term term);
    void addConstraint(ConstraintDescriptor constraint);
    void addClosure(std::string closure);
    void addDependency(std::string dependency);
    void requestTransformation(TransformationDescriptor descriptor);
    /// @brief 该 transformation 是否已被 composition 请求（用于避免重复请求）。
    bool requestsTransformation(std::string_view id) const;
    void addExecutionPolicy(LegacyExecutionPolicy policy);
    void extendExecutionPolicy(
        const std::string& policyId,
        const std::string& equation,
        const std::string& unknown);
    void applyModification(SystemModification modification);
    void applyContribution(SystemContribution contribution);

    void addExecution(ExecutionScope value) { execution.push_back(std::move(value)); }
    void bindNumerics(NumericalBinding value) { numerics.push_back(std::move(value)); }
    std::vector<ExecutionScope> execution;
    std::vector<ExecutionScope> legacyExecution;
    std::vector<NumericalBinding> numerics;

    const RawEquationSystem& rawSystem() const { return system_; }

private:
    RawEquationSystem& system_;
    std::vector<TransformationDescriptor>& transformations_;
    std::vector<LegacyExecutionPolicy>& policies_;
    Provenance origin_;
};

} // namespace SF::System
