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
        std::vector<ExecutionPolicy>& policies,
        Provenance origin);

    void recordContribution(std::string id, std::string name);
    void addUnknown(UnknownDescriptor unknown);
    void addEquation(
        EquationDescriptor descriptor, Equation::Definition definition);
    void extendEquation(const std::string& equationId, Equation::Term term);
    void addConstraint(ConstraintDescriptor constraint);
    void addClosure(std::string closure);
    void addDependency(std::string dependency);
    void requestTransformation(TransformationDescriptor descriptor);
    /// @brief 该 transformation 是否已被 composition 请求（用于避免重复请求）。
    bool requestsTransformation(std::string_view id) const;
    void addExecutionPolicy(ExecutionPolicy policy);
    void extendExecutionPolicy(
        const std::string& policyId,
        const std::string& equation,
        const std::string& unknown);
    void applyModification(SystemModification modification);
    void applyContribution(SystemContribution contribution);

    const RawEquationSystem& rawSystem() const { return system_; }

private:
    RawEquationSystem& system_;
    std::vector<TransformationDescriptor>& transformations_;
    std::vector<ExecutionPolicy>& policies_;
    Provenance origin_;
};

} // namespace SF::System
