#pragma once

/// @file SF_systemContribution.h
/// @brief Model 对方程/约束/闭合/执行策略的中立值贡献。

#include "core/system/SF_equationIR.h"
#include "core/system/SF_solveProgram.h"
#include "core/system/SF_transformationTypes.h"

#include <utility>

namespace SF::System {

struct ContributedEquation {
    EquationDescriptor descriptor;
    Equation::Definition definition;
};

struct EquationExtension {
    std::string equation;
    Equation::Term term;
};

struct PolicyExtension {
    std::string policy;
    std::string equation;
    std::string unknown;
};

/// @brief Model 只产生值；solver/system compiler 负责验证并并入 raw system。
class SystemContribution {
public:
    explicit SystemContribution(const RawEquationSystem& raw) : raw_(raw) {}

    const RawEquationSystem& rawSystem() const { return raw_; }
    void recordContribution(std::string id, std::string name) {
        records.push_back({std::move(id),std::move(name),{}});
    }
    void addUnknown(UnknownDescriptor value) { unknowns.push_back(std::move(value)); }
    void addEquation(EquationDescriptor descriptor, Equation::Definition definition) {
        equations.push_back({std::move(descriptor),std::move(definition)});
    }
    void extendEquation(std::string id, Equation::Term term) {
        equationExtensions.push_back({std::move(id),std::move(term)});
    }
    void addConstraint(ConstraintDescriptor value) { constraints.push_back(std::move(value)); }
    void addClosure(std::string value) { closures.push_back(std::move(value)); }
    void addDependency(std::string value) { dependencies.push_back(std::move(value)); }
    void requestTransformation(TransformationDescriptor value) {
        transformations.push_back(std::move(value));
    }
    void addExecutionPolicy(ExecutionPolicy value) { policies.push_back(std::move(value)); }
    void extendExecutionPolicy(std::string policy, std::string equation,
                               std::string unknown) {
        policyExtensions.push_back({std::move(policy),std::move(equation),
                                    std::move(unknown)});
    }

    std::vector<ContributionRecord> records;
    std::vector<UnknownDescriptor> unknowns;
    std::vector<ContributedEquation> equations;
    std::vector<EquationExtension> equationExtensions;
    std::vector<ConstraintDescriptor> constraints;
    std::vector<std::string> closures;
    std::vector<std::string> dependencies;
    std::vector<TransformationDescriptor> transformations;
    std::vector<ExecutionPolicy> policies;
    std::vector<PolicyExtension> policyExtensions;

private:
    const RawEquationSystem& raw_;
};

} // namespace SF::System
