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
        std::vector<LegacyExecutionPolicy>& policies,
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

void SystemCompositionBuilder::addState(StateSymbol unknown) {
    if (unknown.id.empty()) throw std::runtime_error("Unknown contribution requires an id.");
    if (hasUnknown(system_,unknown.id)) {
        throw std::runtime_error("Duplicate unknown contribution '"+unknown.id+"'.");
    }
    unknown.origin = origin_;
    system_.state.add(std::move(unknown));
}

void SystemCompositionBuilder::addEquation(
        EquationDescriptor descriptor, SF::Equation::Definition definition) {
    if (descriptor.id != definition.name) {
        throw std::runtime_error(
            "Equation id does not match definition '"+descriptor.id+"'.");
    }
    if (hasEquation(system_,descriptor.id)) {
        throw std::runtime_error("Duplicate equation contribution '"+descriptor.id+"'.");
    }
    descriptor.origin = origin_;
    system_.registry.add(formulaFromEquation(definition,descriptor.origin));
    system_.legacyEquations.push_back(std::move(descriptor));
    system_.legacyDefinitions.add(std::move(definition));
}

void SystemCompositionBuilder::addEquation(
        EquationDescriptor descriptor, SF::Equation::Definition definition,
        Equation formula) {
    if (descriptor.id!=definition.name || descriptor.id!=formula.id)
        throw std::runtime_error("Authored Equation and equation IDs differ.");
    if (hasEquation(system_,descriptor.id))
        throw std::runtime_error("Duplicate equation contribution '"+descriptor.id+"'.");
    descriptor.origin=origin_;
    formula.origin=origin_;
    formula.authored=true;
    system_.registry.add(std::move(formula));
    system_.legacyEquations.push_back(std::move(descriptor));
    system_.legacyDefinitions.add(std::move(definition));
}

void SystemCompositionBuilder::addEquation(Equation equation) {
    equation.origin=origin_;
    equation.authored=true;
    system_.registry.add(std::move(equation));
}

void SystemCompositionBuilder::extendEquation(
        const std::string& equationId, SF::Equation::Term term) {
    if (!hasEquation(system_,equationId)) {
        throw std::runtime_error(
            "Cannot extend undeclared equation '"+equationId+"'.");
    }
    auto equation=system_.registry.at(equationId);
    const auto legacy=std::find_if(system_.legacyEquations.begin(),system_.legacyEquations.end(),
        [&](const EquationDescriptor& value) { return value.id==equationId; });
    if (equation.authored) {
        if (term.kind!=SF::Equation::TermKind::Source)
            throw std::runtime_error("Native equation extension requires an explicit mathematical AST contribution.");
        equation.rhs=FormulaExpr::add(std::move(equation.rhs),
            FormulaExpr::op("source",{FormulaExpr::symbol(term.primary.name)},term.primary.name));
        system_.registry.replace(std::move(equation));
        if (legacy!=system_.legacyEquations.end())
            system_.legacyDefinitions.addRightTerm(equationId,std::move(term));
    } else {
        system_.legacyDefinitions.addRightTerm(equationId,std::move(term));
        system_.registry.replace(formulaFromEquation(
            system_.legacyDefinitions.at(equationId),origin_));
    }
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

void SystemCompositionBuilder::addExecutionPolicy(LegacyExecutionPolicy policy) {
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
        [&](const LegacyExecutionPolicy& item) { return item.id == policyId; });
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
    if (contribution.immersed) {
        if (system_.immersed) throw std::runtime_error("Duplicate resolved immersed port selection.");
        system_.immersed=contribution.immersed;
    }
    placement.insert(placement.end(),contribution.placement.begin(),contribution.placement.end());
    for (auto& node:contribution.execution) {
        node.origin=origin_; execution.push_back(std::move(node));
    }
    for (auto& node:contribution.legacyExecution) {
        node.origin=origin_; legacyExecution.push_back(std::move(node));
    }
    for (auto& binding:contribution.numerics) numerics.push_back(std::move(binding));
    for (auto& unknown : contribution.states) addState(std::move(unknown));
    for (auto& id:contribution.requiredStates) requireState(std::move(id));
    for (auto& equation : contribution.legacyEquations) {
        addEquation(std::move(equation.descriptor),std::move(equation.definition));
    }
    for (auto& formula:contribution.registeredEquations) addEquation(std::move(formula));
    for (const auto& extension:contribution.mathematicalExtensions) {
        bool matched=false;
        auto entries=system_.registry.entries();
        for (auto equation:entries) {
            const bool applies=std::any_of(extension.families.begin(),extension.families.end(),
                [&](const std::string& family) {
                    return equation.id==family || equation.id.rfind(family+".",0)==0;
                });
            if (!applies) continue;
            equation.rhs=FormulaExpr::add(std::move(equation.rhs),extension.rhs);
            equation.authored=true;
            system_.registry.replace(std::move(equation));
            matched=true;
        }
        if (!matched) throw std::runtime_error("Source contribution has no registered target equation.");
    }
    for (auto& extension : contribution.equationExtensions) {
        if (extension.role) {
            bool matched=false;
            for (const auto& equation : system_.legacyEquations)
                if (equation.role == *extension.role
                    || (extension.alternativeRole
                        && equation.role==*extension.alternativeRole)) {
                    extendEquation(equation.id,extension.term);
                    matched=true;
                }
            if (!matched)
                throw std::runtime_error(
                    "Model term '"+extension.term.primary.name
                    +"' has no applicable equation role in the composed system.");
        } else {
            extendEquation(extension.equation,std::move(extension.term));
        }
    }
    for (auto& constraint : contribution.constraints) addConstraint(std::move(constraint));
    for (auto& closure : contribution.boundaryClosures) {
        if (closure.id.empty() || closure.provider.empty() || closure.reads.empty()
            || closure.writes.empty() || closure.order.empty())
            throw std::runtime_error("Incomplete boundary closure contract.");
        if (std::any_of(system_.boundaryClosures.begin(),system_.boundaryClosures.end(),
            [&](const auto& previous) { return previous.id==closure.id; }))
            throw std::runtime_error("Duplicate boundary closure contract: "+closure.id);
        closure.origin=origin_;
        system_.boundaryClosures.push_back(std::move(closure));
    }
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
