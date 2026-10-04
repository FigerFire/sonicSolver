/// @file SF_systemValidator.cpp
/// @brief 启动阶段统一拒绝未声明未知量、悬空约束和缺失执行能力。

#include "SF_systemValidator.h"
#include "SF_solvePlan.h"

#include <map>
#include <set>
#include <stdexcept>

namespace SF::System {
namespace {

void validateComposition(const RawEquationSystem& raw) {
    std::set<std::string> contributions;
    for (const auto& contribution : raw.contributions) {
        if (contribution.id.empty()
            || !contributions.insert(contribution.id).second) {
            throw std::runtime_error(
                "Composition has an invalid/duplicate contribution.");
        }
    }
    std::set<std::string> unknowns;
    for (const auto& unknown : raw.state.symbols()) {
        if (unknown.id.empty() || !unknowns.insert(unknown.id).second) {
            throw std::runtime_error(
                "Raw equation system has an invalid/duplicate unknown.");
        }
        if (unknown.id == "pPrime") {
            throw std::runtime_error(
                "Raw equation system contains algorithmic pressure-correction "
                "workspace 'pPrime'.");
        }
    }
    std::set<std::string> equations;
    for (const auto& equation : raw.legacyEquations) {
        if (equation.id.empty() || !equations.insert(equation.id).second
            || raw.legacyDefinitions.at(equation.id).name != equation.id) {
            throw std::runtime_error(
                "Raw equation system has an invalid/duplicate equation.");
        }
        if (equation.category == EquationCategory::AlgorithmicDerivedEquation) {
            throw std::runtime_error(
                "Raw equation system contains algorithmic-derived equation '"
                +equation.id+"'.");
        }
    }
}

void validateCompiledEquations(const ExecutableEquationSystem& executable) {
    std::set<std::string> ids;
    for (const auto& compiled : executable.compiledEquations) {
        if (compiled.equationId.empty()
            || !ids.insert(compiled.equationId).second) {
            throw std::runtime_error(
                "Executable system has an invalid/duplicate compiled equation.");
        }
        if (!executable.registry.contains(compiled.equationId))
            throw std::runtime_error("Compiled equation has no registry definition: "+compiled.equationId);
        if (compiled.origin.kind != OriginKind::Generated
            || compiled.resources.empty()) {
            throw std::runtime_error(
                "Compiled equation '"+compiled.equationId
                +"' has no generated storage/operator binding.");
        }
        for (const auto& resource : compiled.resources) {
            if (resource.symbol.empty() || resource.storage.empty()
                || resource.components <= 0 || resource.componentOffset < 0) {
                throw std::runtime_error(
                    "Compiled equation '"+compiled.equationId
                    +"' has an invalid resource binding.");
            }
        }
    }
}

void validateTransformations(const ResolvedSimulationSystem& system) {
    std::set<std::string> ids;
    for (const auto& transformation : system.transformations) {
        if (transformation.descriptor.id.empty()
            || !ids.insert(transformation.descriptor.id).second) {
            throw std::runtime_error(
                "Transformation pipeline has an invalid/duplicate request.");
        }
        if (transformation.state == TransformationState::NotRegistered
            || transformation.state == TransformationState::RegisteredButInvalid
            || transformation.state == TransformationState::RegisteredAndActive) {
            throw std::runtime_error(
                "Transformation '"+transformation.descriptor.id
                +"' did not reach a terminal valid state.");
        }
        if (transformation.descriptor.explicitlyRequested
            && transformation.state != TransformationState::Applied) {
            throw std::runtime_error(
                "Explicit transformation '"+transformation.descriptor.id
                +"' was not applied: "+transformation.reason+".");
        }
    }
}

void validatePlanNode(const SolvePlanNode& node) {
    if (node.id.empty()) {
        throw std::runtime_error("Compiled solve plan contains an unnamed node.");
    }
    if (node.children.empty() && node.operation.empty()) {
        throw std::runtime_error(
            "Compiled solve-plan leaf '"+node.id
            +"' has no explicit runtime operation ID.");
    }
    for (const auto& child : node.children) validatePlanNode(child);
}

} // namespace

void validate(const ResolvedSimulationSystem& system) {
    if (system.coupling.status==CouplingStatus::Invalid)
        throw std::runtime_error("Invalid explicitly selected coupling preset: "+system.coupling.reason);
    validateComposition(system.rawSystem);
    validateTransformations(system);
    validatePlanNode(system.solvePlan.root);
    const auto& executable = system.executableSystem;
    validateCompiledEquations(executable);
    std::set<std::string> providers;
    for (const auto& provider : system.runtime.providerRequirements) {
        if (provider.id.empty() || provider.responsibility.empty()
            || !providers.insert(provider.id).second) {
            throw std::runtime_error(
                "Resolved system has an invalid/duplicate provider requirement.");
        }
    }
    std::set<std::string> runtimeServices;
    for (const auto& service : system.runtime.runtimeServiceRequirements) {
        if (service.id.empty() || service.reason.empty()
            || !runtimeServices.insert(service.id).second) {
            throw std::runtime_error(
                "Resolved system has an invalid/duplicate runtime-service requirement.");
        }
    }
    const bool hasEulerianProvider =
        providers.count("flow.eulerian-pressure") != 0;
    const bool hasConservativeProvider =
        providers.count("flow.conservative") != 0;
    const bool hasPressureProvider=providers.count("flow.pressure-operators")!=0;
    if ((int)hasEulerianProvider+(int)hasConservativeProvider+(int)hasPressureProvider>1
        || (providers.count("flow.rhie-chow") && !hasPressureProvider)) {
        throw std::runtime_error(
            "Resolved system cannot bind two flow operation provider groups.");
    }
    std::set<OpId> resolvedBindings;
    std::set<OpId> missingBindings;
    for (const auto& binding : system.runtime.operationBindings) {
        if (binding.operation.empty()
            || !resolvedBindings.insert(binding.operation).second) {
            throw std::runtime_error(
                "Provider resolution has an invalid/duplicate operation binding.");
        }
        if (binding.status == BindingStatus::Resolved) {
            if (binding.provider.empty() || !providers.count(binding.provider)) {
                throw std::runtime_error(
                    "Resolved operation provider is absent from runtime requirements.");
            }
        } else {
            if (!binding.provider.empty() || binding.reason.empty()) {
                throw std::runtime_error(
                    "Unsupported operation binding lacks an exact reason.");
            }
            missingBindings.insert(binding.operation);
        }
    }
    const auto requiredOperations = SolvePlanner::requiredOperations(system.solvePlan);
    if (requiredOperations != system.runtime.report.requiredOperations
        || resolvedBindings.size() != requiredOperations.size()) {
        throw std::runtime_error(
            "Operation bindings do not cover the compiled solve plan.");
    }
    for (const OpId& id : requiredOperations) {
        if (!resolvedBindings.count(id)) {
            throw std::runtime_error(
                "Plan operation '"+id+"' has no provider resolution result.");
        }
    }
    const auto validateOwners=[&](const auto& self,const SolvePlanNode& node)->void {
        if (!node.operation.empty()) {
            const auto binding=std::find_if(system.runtime.operationBindings.begin(),system.runtime.operationBindings.end(),
                [&](const auto& item) { return item.operation==node.operation; });
            if (binding!=system.runtime.operationBindings.end() && binding->status==BindingStatus::Resolved
                && (node.provider.empty() || node.provider!=binding->provider))
                throw std::runtime_error("Compiled operation ownership disagrees with its runtime binding: "+node.operation);
        }
        for (const auto& child:node.children) self(self,child);
    };
    validateOwners(validateOwners,system.solvePlan.root);
    if (std::set<OpId>(system.runtime.report.missingOperations.begin(),
                       system.runtime.report.missingOperations.end())
            != missingBindings
        || system.runtime.report.missingOperations.size()
            != missingBindings.size()) {
        throw std::runtime_error(
            "RuntimeReport missing operations disagree with provider bindings.");
    }
    std::set<std::string> unknowns;
    for (const auto& unknown : executable.state.symbols()) {
        if (unknown.id.empty() || unknown.name.empty()
            || unknown.components <= 0
            || (unknown.storageBinding != StorageBinding::SpecializedExecutor
                && unknown.runtimeStorageRequired
                && unknown.storageKey.empty())
            || !unknowns.insert(unknown.id).second) {
            throw std::runtime_error(
                "Resolved system has an invalid/duplicate unknown.");
        }
    }
    std::set<std::string> equations;
    for (const auto& equation:executable.registry.entries()) equations.insert(equation.id);
    std::vector<std::string> equationIds;
    for (const auto& equation : executable.legacyEquations) {
        if (equation.id.empty() || equation.name.empty()
            || !executable.registry.contains(equation.id)) {
            throw std::runtime_error(
                "Resolved system has an invalid/duplicate equation.");
        }
        const auto& definition = equationDefinition(system,equation.id);
        equationIds.push_back(equation.id);
        if (definition.name != equation.id) {
            throw std::runtime_error(
                "Equation descriptor and executable definition disagree for '"
                + equation.id + "'.");
        }
        for (const std::string& unknown : equation.solvedUnknowns) {
            if (unknowns.find(unknown) == unknowns.end()) {
                throw std::runtime_error(
                    "Equation '"+equation.id
                    +"' references undeclared unknown '"+unknown+"'.");
            }
        }
    }
    std::set<std::string> constraints;
    for (const auto& constraint : executable.constraints) {
        if (constraint.id.empty() || constraint.equation.empty()
            || (!constraint.multiplierUnknown.empty()
                && unknowns.find(constraint.multiplierUnknown) == unknowns.end())
            || !constraints.insert(constraint.id).second) {
            throw std::runtime_error(
                "Resolved system has an invalid multiplier constraint.");
        }
    }
    std::set<std::string> blocks;
    std::map<std::string,int> equationBlockReferences;
    for (const auto& block : system.solvePlan.blocks) {
        if (block.id.empty() || block.strategy.empty()
            || !blocks.insert(block.id).second) {
            throw std::runtime_error(
                "Resolved system has an invalid/duplicate solve block.");
        }
        if (block.strategyKind == FDM::SolveStrategyKind::AlgebraicUpdate) {
            throw std::runtime_error(
                "Solve block '"+block.id
                +"' has no executable typed solve strategy.");
        }
        for (const auto& equation : block.equations) {
            if (equations.find(equation) == equations.end()) {
                throw std::runtime_error(
                    "Solve block '"+block.id
                    +"' references undeclared equation '"+equation+"'.");
            }
            ++equationBlockReferences[equation];
        }
        for (const auto& constraint : block.constraints) {
            if (constraints.find(constraint) == constraints.end()) {
                throw std::runtime_error(
                    "Solve block '"+block.id
                    +"' references undeclared constraint '"+constraint+"'.");
            }
        }
        for (const auto& unknown : block.unknowns) {
            if (unknowns.find(unknown) == unknowns.end()) {
                throw std::runtime_error(
                    "Solve block '"+block.id
                    +"' references undeclared unknown '"+unknown+"'.");
            }
        }
    }
    // Unmigrated compatibility equations still require explicit execution ownership.
    // Native occurrences were already validated by their selected provider.
    if (!equationIds.empty() && system.solvePlan.compiledProgram.steps.empty()) {
        for (const auto& equation:equationIds) {
            (void)executable.legacyDefinitions.at(equation);
            if (equationBlockReferences[equation]==0 && system.runtime.report.status==RuntimeStatus::Runnable)
                throw std::runtime_error("Legacy equation lacks execution ownership: "+equation);
        }
    }
    for (const auto& call:system.solvePlan.compiledProgram.steps) {
        (void)executable.registry.at(call.source.equation);
        if (call.target.resources.empty())
            throw std::runtime_error("Compiled occurrence lacks a typed target binding.");
    }
    std::set<std::string> capabilities;
    for (const auto& requirement : system.runtime.requirements) {
        // required / available / reason 必须语义一致：name 声明能力，
        // available 只表示"本 binary 提供该能力"，detail 说明它是什么。
        if (requirement.name.empty() || requirement.detail.empty()) {
            throw std::runtime_error(
                "Resolved system has an unnamed or unjustified runtime "
                "capability requirement.");
        }
        if (!capabilities.insert(requirement.name).second) {
            throw std::runtime_error(
                "Resolved system repeats capability '"+requirement.name+"'.");
        }
        if (requirement.required && !requirement.available) {
            throw std::runtime_error(
                "Resolved system requires unavailable capability '"
                +requirement.name+"': "+requirement.detail+".");
        }
        if (!requirement.available
            && system.runtime.report.status == RuntimeStatus::Runnable) {
            // 不允许 "input says A / compiled says B / runtime silently runs C"：
            // 任何 unavailable capability 都必须让系统不可执行。
            throw std::runtime_error(
                "Resolved system is reported Runnable while capability '"
                +requirement.name+"' is unavailable: "+requirement.detail+".");
        }
    }
    std::set<std::string> workspaces;
    for (const auto& workspace : system.runtime.workspaceRequirements) {
        if (workspace.id.empty() || workspace.components <= 0
            || !workspaces.insert(workspace.id).second) {
            throw std::runtime_error(
                "Resolved system has an invalid/duplicate workspace requirement.");
        }
    }
}

} // namespace SF::System
