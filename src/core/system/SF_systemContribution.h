#pragma once

/// @file SF_systemContribution.h
/// @brief Model 对方程/约束/闭合/执行策略的中立值贡献。

#include "core/system/SF_equationIR.h"
#include "core/system/SF_solveProgram.h"
#include "core/system/SF_transformationTypes.h"

#include <optional>
#include <utility>

namespace SF::System {

struct ContributedEquation {
    EquationDescriptor descriptor;
    SF::Equation::Definition definition;
};

struct EquationExtension {
    std::string equation;
    std::optional<LegacyEquationRole> role;
    std::optional<LegacyEquationRole> alternativeRole;
    SF::Equation::Term term;
};

/// @brief Extend one mathematical family, including explicitly namespaced phase equations.
struct MathematicalExtension {
    std::vector<std::string> families;
    FormulaExpr rhs;
};

struct PolicyExtension {
    std::string policy;
    std::string equation;
    std::string unknown;
};

struct ContributionProviderRequirement {
    std::string id;
    std::string responsibility;
};

/// @brief Model 只产生值；solver/system compiler 负责验证并并入 raw system。
class SystemContribution {
public:
    SystemContribution() = default;
    void recordContribution(std::string id, std::string name) {
        records.push_back({std::move(id),std::move(name),{}});
    }
    void addState(StateSymbol value) { states.push_back(std::move(value)); }
    void requireState(std::string id) { requiredStates.push_back(std::move(id)); }
    void addEquation(EquationDescriptor descriptor, SF::Equation::Definition definition) {
        legacyEquations.push_back({std::move(descriptor),std::move(definition)});
    }
    /// @brief Add a mathematical formula without requiring a legacy
    /// SF::Equation::Definition or semantic LegacyEquationRole.
    void addEquation(Equation value) { registeredEquations.push_back(std::move(value)); }
    void extendMathematics(std::vector<std::string> families,FormulaExpr rhs) {
        mathematicalExtensions.push_back({std::move(families),std::move(rhs)});
    }
    void extendEquation(std::string id, SF::Equation::Term term) {
        equationExtensions.push_back(
            {std::move(id),std::nullopt,std::nullopt,std::move(term)});
    }
    /// @brief Extend every equation with this mathematical role, independent of preset IDs.
    void extendEquationRole(LegacyEquationRole role, SF::Equation::Term term) {
        equationExtensions.push_back({{},role,std::nullopt,std::move(term)});
    }
    void extendEquationEitherRole(
            LegacyEquationRole first,LegacyEquationRole second,SF::Equation::Term term) {
        equationExtensions.push_back({{},first,second,std::move(term)});
    }
    void addConstraint(ConstraintDescriptor value) { constraints.push_back(std::move(value)); }
    void addClosure(std::string value) { closures.push_back(std::move(value)); }
    void requireProvider(std::string id,std::string responsibility) {
        providerRequirements.push_back({std::move(id),std::move(responsibility)});
    }
    void addDependency(std::string value) { dependencies.push_back(std::move(value)); }
    void requestTransformation(TransformationDescriptor value) {
        transformations.push_back(std::move(value));
    }
    void addExecutionPolicy(LegacyExecutionPolicy value) { policies.push_back(std::move(value)); }
    void extendExecutionPolicy(std::string policy, std::string equation,
                               std::string unknown) {
        policyExtensions.push_back({std::move(policy),std::move(equation),
                                    std::move(unknown)});
    }

    /// @brief Modules contribute scope-local calls and numerical defaults independently.
    void addExecution(ExecutionScope value) { execution.push_back(std::move(value)); }
    void bindNumerics(NumericalBinding value) { numerics.push_back(std::move(value)); }
    /// @brief Migration boundary: declarations whose legacy backend still owns local topology.
    void addLegacyExecution(ExecutionScope value) { legacyExecution.push_back(std::move(value)); }
    std::vector<ExecutionScope> execution;
    std::vector<ExecutionScope> legacyExecution;
    std::vector<NumericalBinding> numerics;

    std::vector<ContributionRecord> records;
    std::vector<StateSymbol> states;
    std::vector<std::string> requiredStates;
    std::vector<ContributedEquation> legacyEquations;
    std::vector<Equation> registeredEquations;
    std::vector<EquationExtension> equationExtensions;
    std::vector<MathematicalExtension> mathematicalExtensions;
    std::vector<ConstraintDescriptor> constraints;
    std::vector<std::string> closures;
    std::vector<ContributionProviderRequirement> providerRequirements;
    std::vector<std::string> dependencies;
    std::vector<TransformationDescriptor> transformations;
    std::vector<LegacyExecutionPolicy> policies;
    std::vector<PolicyExtension> policyExtensions;

};

} // namespace SF::System
