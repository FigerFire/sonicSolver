#pragma once

/// @file SF_equationContribution.h
/// @brief Builtin、model 与 user 共用的 equation-system composition contract。

#include "SF_resolvedSimulationSystem.h"
#include "core/system/SF_systemContribution.h"
#include "SF_pressureCoupling.h"

#include <memory>
#include <vector>

namespace SF::System {

/// @brief 向 raw system 注册普通数学对象，不执行 transformation 或 timestep。
class SystemCompositionBuilder {
public:
    SystemCompositionBuilder(
        RawEquationSystem& system,
        std::vector<TransformationDescriptor>& transformations,
        Provenance origin);

    void recordContribution(std::string id, std::string name,std::string module={});
    void addState(StateSymbol unknown);
    /// A request activates catalog metadata at freeze; it adds no equations/storage.
    void requireState(std::string id) { requiredStates.push_back(std::move(id)); }
    std::vector<std::string> requiredStates;
    void addEquation(Equation equation);
    void addConstraint(ConstraintDescriptor constraint);
    void addClosure(std::string closure,std::vector<std::string> equations={},std::vector<std::string> states={},std::vector<std::string> outputs={});
    void addDependency(std::string dependency);
    void requestTransformation(TransformationDescriptor descriptor);
    /// @brief 该 transformation 是否已被 composition 请求（用于避免重复请求）。
    bool requestsTransformation(std::string_view id) const;
    void applyModification(SystemModification modification);
    void applyContribution(SystemContribution contribution);
    /// @brief 所有模块/用户注册完成后合并显式加性项，不推断系数组合或重排浮点项。
    void finalizeMathematics();

    void addExecution(ExecutionScope value) { execution.push_back(std::move(value)); }
    void bindNumerics(NumericalBinding value) { value.origin=origin_;numerics.push_back(std::move(value)); }
    std::optional<PressureCouplingInputs> pressureInputs;
    std::vector<PlacementRequirement> placement;
    std::vector<ExecutionScope> execution;
    std::vector<NumericalBinding> numerics;

    const RawEquationSystem& rawSystem() const { return system_; }

private:
    RawEquationSystem& system_;
    std::vector<TransformationDescriptor>& transformations_;
    Provenance origin_;
    std::vector<std::pair<MathematicalExtension,Provenance>> extensions_;
};

} // namespace SF::System
