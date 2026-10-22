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
        Provenance origin)
    : system_(system),
      transformations_(transformations),
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
    if(unknown.origin.source.empty())unknown.origin = origin_;
    system_.state.add(std::move(unknown));
}

void SystemCompositionBuilder::addEquation(Equation equation) {
    system_.state.validateSymbols(equation.lhs);system_.state.validateSymbols(equation.rhs);
    if(equation.origin.source.empty())equation.origin=origin_;
    system_.registry.add(std::move(equation));
}

void SystemCompositionBuilder::addConstraint(ConstraintDescriptor constraint) {
    if (hasConstraint(system_,constraint.id)) {
        throw std::runtime_error(
            "Duplicate constraint contribution '"+constraint.id+"'.");
    }
    constraint.origin = origin_;
    system_.constraints.push_back(std::move(constraint));
}

void SystemCompositionBuilder::addClosure(std::string closure,
        std::vector<std::string> equations,std::vector<std::string> states) {
    if (closure.empty()) throw std::runtime_error("Closure contribution is empty.");
    if (std::any_of(system_.closures.begin(),system_.closures.end(),[&](const auto& item) { return item.description==closure; }))
        throw std::runtime_error("Duplicate closure contribution: "+closure);
    system_.closures.push_back({std::move(closure),std::move(equations),std::move(states),origin_});
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
    const auto parentOrigin=origin_;
    if (!contribution.records.empty()) {
        origin_.source+=" /";
        for (const auto& record:contribution.records) origin_.source+=" "+record.id;
    }
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
    for (auto& binding:contribution.numerics) numerics.push_back(std::move(binding));
    for (auto& unknown : contribution.states) addState(std::move(unknown));
    for (auto& id:contribution.requiredStates) requireState(std::move(id));
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
            const auto identity=canonicalFormula(Equation{{},FormulaExpr::constantValue(0),extension.rhs,{}});
            if (std::any_of(system_.modifications.begin(),system_.modifications.end(),[&](const auto& item) {
                    return item.kind==ModificationKind::Extend && item.targetKind=="equation"
                        && item.targetId==equation.id && item.detail==identity;
                })) throw std::runtime_error("Duplicate mathematical contribution to equation: "+equation.id);
            system_.modifications.push_back({ModificationKind::Extend,"equation",equation.id,identity,origin_});
            equation.rhs=FormulaExpr::add(std::move(equation.rhs),extension.rhs);
            system_.registry.replace(std::move(equation));
            matched=true;
        }
        if (!matched) throw std::runtime_error("Source contribution has no registered target equation.");
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
    for (auto& closure : contribution.closures)
        addClosure(std::move(closure.description),std::move(closure.equations),std::move(closure.states));
    for (auto& dependency : contribution.dependencies) addDependency(std::move(dependency));
    for (auto& descriptor : contribution.transformations) {
        requestTransformation(std::move(descriptor));
    }
    origin_=parentOrigin;
}

} // namespace SF::System
