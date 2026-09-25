/// @file SF_equationContribution.cpp
/// @brief Unified system composition 与显式 modification 语义。

#include "SF_equationContribution.h"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace SF::System {

SystemCompositionBuilder::SystemCompositionBuilder(
        RawEquationSystem& system,
        std::vector<TransformationDescriptor>& transformations,
        std::vector<ExecutionPolicy>& policies,
        Provenance origin)
    : system_(system),
      transformations_(transformations),
      policies_(policies),
      origin_(std::move(origin)) {}

void SystemCompositionBuilder::recordContribution(
        std::string id, std::string name) {
    if (id.empty()) throw std::runtime_error("System contribution requires an id.");
    const auto duplicate = std::find_if(
        system_.contributions.begin(),system_.contributions.end(),
        [&](const ContributionRecord& item) { return item.id == id; });
    if (duplicate != system_.contributions.end()) {
        throw std::runtime_error("Duplicate system contribution '"+id+"'.");
    }
    system_.contributions.push_back({std::move(id),std::move(name),origin_});
}

void SystemCompositionBuilder::addUnknown(UnknownDescriptor unknown) {
    if (unknown.id.empty()) throw std::runtime_error("Unknown contribution requires an id.");
    if (hasUnknown(system_,unknown.id)) {
        throw std::runtime_error("Duplicate unknown contribution '"+unknown.id+"'.");
    }
    unknown.origin = origin_;
    system_.unknowns.push_back(std::move(unknown));
}

void SystemCompositionBuilder::addEquation(
        EquationDescriptor descriptor, Equation::Definition definition) {
    if (descriptor.id != definition.name) {
        throw std::runtime_error(
            "Equation id does not match definition '"+descriptor.id+"'.");
    }
    if (hasEquation(system_,descriptor.id)) {
        throw std::runtime_error("Duplicate equation contribution '"+descriptor.id+"'.");
    }
    descriptor.origin = origin_;
    system_.equations.push_back(std::move(descriptor));
    system_.equationDefinitions.add(std::move(definition));
}

void SystemCompositionBuilder::extendEquation(
        const std::string& equationId, Equation::Term term) {
    if (!hasEquation(system_,equationId)) {
        throw std::runtime_error(
            "Cannot extend undeclared equation '"+equationId+"'.");
    }
    system_.equationDefinitions.addRightTerm(equationId,std::move(term));
    system_.modifications.push_back({
        ModificationKind::Extend,"equation",equationId,
        "append right-hand term",origin_});
}

void SystemCompositionBuilder::addConstraint(ConstraintDescriptor constraint) {
    if (hasConstraint(system_,constraint.id)) {
        throw std::runtime_error(
            "Duplicate constraint contribution '"+constraint.id+"'.");
    }
    constraint.origin = origin_;
    system_.constraints.push_back(std::move(constraint));
}

void SystemCompositionBuilder::addClosure(std::string closure) {
    if (closure.empty()) throw std::runtime_error("Closure contribution is empty.");
    system_.closures.push_back(std::move(closure));
}

void SystemCompositionBuilder::addDependency(std::string dependency) {
    if (dependency.empty()) throw std::runtime_error("Dependency contribution is empty.");
    system_.dependencies.push_back(std::move(dependency));
}

bool SystemCompositionBuilder::requestsTransformation(
        std::string_view id) const {
    return std::any_of(
        transformations_.begin(),transformations_.end(),
        [&](const TransformationDescriptor& descriptor) {
            return descriptor.id == id;
        });
}

void SystemCompositionBuilder::requestTransformation(
        TransformationDescriptor descriptor) {
    if (descriptor.id.empty()) {
        throw std::runtime_error("Transformation request requires an id.");
    }
    descriptor.origin = origin_;
    transformations_.push_back(std::move(descriptor));
}

void SystemCompositionBuilder::addExecutionPolicy(ExecutionPolicy policy) {
    if (policy.id.empty()) throw std::runtime_error("Execution policy requires an id.");
    policy.origin = origin_;
    policies_.push_back(std::move(policy));
}

void SystemCompositionBuilder::extendExecutionPolicy(
        const std::string& policyId,
        const std::string& equation,
        const std::string& unknown) {
    const auto found = std::find_if(
        policies_.begin(),policies_.end(),
        [&](const ExecutionPolicy& item) { return item.id == policyId; });
    if (found == policies_.end()) {
        throw std::runtime_error(
            "Cannot extend undeclared execution policy '"+policyId+"'.");
    }
    found->equations.push_back(equation);
    found->unknowns.push_back(unknown);
}

void SystemCompositionBuilder::applyModification(SystemModification modification) {
    modification.origin = origin_;
    if (modification.kind == ModificationKind::Add) {
        throw std::runtime_error(
            "SystemModification::Add must use a typed add operation.");
    }
    if (modification.kind == ModificationKind::Extend) {
        throw std::runtime_error(
            "SystemModification::Extend must provide a typed equation term.");
    }
    // Replace/disable 需要完整 payload 或类型化 target API。拒绝不完整请求，
    // 避免重新引入同名后写覆盖前写的隐式规则。
    throw std::runtime_error(
        "System modification '"+std::string(toString(modification.kind))
        +"' for "+modification.targetKind+" '"+modification.targetId
        +"' is registered but lowering is not implemented.");
}

void SystemCompositionBuilder::applyContribution(SystemContribution contribution) {
    for (auto& record : contribution.records) {
        recordContribution(std::move(record.id),std::move(record.name));
    }
    for (auto& unknown : contribution.unknowns) addUnknown(std::move(unknown));
    for (auto& equation : contribution.equations) {
        addEquation(std::move(equation.descriptor),std::move(equation.definition));
    }
    for (auto& extension : contribution.equationExtensions) {
        extendEquation(extension.equation,std::move(extension.term));
    }
    for (auto& constraint : contribution.constraints) addConstraint(std::move(constraint));
    for (auto& closure : contribution.closures) addClosure(std::move(closure));
    for (auto& dependency : contribution.dependencies) addDependency(std::move(dependency));
    for (auto& descriptor : contribution.transformations) {
        requestTransformation(std::move(descriptor));
    }
    for (auto& policy : contribution.policies) addExecutionPolicy(std::move(policy));
    for (auto& extension : contribution.policyExtensions) {
        extendExecutionPolicy(extension.policy,extension.equation,
                              extension.unknown);
    }
}

} // namespace SF::System
