#pragma once

/// @file SF_systemContribution.h
/// @brief Model 对方程/约束/闭合/执行策略的中立值贡献。

#include "core/system/SF_equationIR.h"
#include "core/system/SF_solveProgram.h"
#include "core/system/SF_transformationTypes.h"

#include <optional>
#include <utility>

namespace SF::System {

/// @brief Extend one mathematical family, including explicitly namespaced phase equations.
struct MathematicalExtension {
    std::vector<std::string> families;
    FormulaExpr rhs;
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
    /// @brief 注册与用户输入共用的数学 AST。
    void addEquation(Equation value) { registeredEquations.push_back(std::move(value)); }
    void extendMathematics(std::vector<std::string> families,FormulaExpr rhs) {
        mathematicalExtensions.push_back({std::move(families),std::move(rhs)});
    }
    void addConstraint(ConstraintDescriptor value) { constraints.push_back(std::move(value)); }
    void addBoundaryClosure(BoundaryClosureContract value) { boundaryClosures.push_back(std::move(value)); }
    void addClosure(std::string value) { closures.push_back(std::move(value)); }
    void requireProvider(std::string id,std::string responsibility) {
        providerRequirements.push_back({std::move(id),std::move(responsibility)});
    }
    void addDependency(std::string value) { dependencies.push_back(std::move(value)); }
    void requestTransformation(TransformationDescriptor value) {
        transformations.push_back(std::move(value));
    }
    /// @brief Modules contribute scope-local calls and numerical defaults independently.
    void addExecution(ExecutionScope value) { execution.push_back(std::move(value)); }
    void bindNumerics(NumericalBinding value) { numerics.push_back(std::move(value)); }
    /// @brief Migration boundary: declarations whose legacy backend still owns local topology.
    std::optional<FDM::ImmersedAlgorithmDescriptor> immersed;
    std::vector<PlacementRequirement> placement;
    std::vector<ExecutionScope> execution;
    std::vector<NumericalBinding> numerics;

    std::vector<ContributionRecord> records;
    std::vector<StateSymbol> states;
    std::vector<std::string> requiredStates;
    std::vector<Equation> registeredEquations;
    std::vector<MathematicalExtension> mathematicalExtensions;
    std::vector<ConstraintDescriptor> constraints;
    std::vector<BoundaryClosureContract> boundaryClosures;
    std::vector<std::string> closures;
    std::vector<ContributionProviderRequirement> providerRequirements;
    std::vector<std::string> dependencies;
    std::vector<TransformationDescriptor> transformations;

};

} // namespace SF::System
