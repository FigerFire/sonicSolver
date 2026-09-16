/// @file SF_systemValidator.cpp
/// @brief 启动阶段统一拒绝未声明未知量、悬空约束和缺失执行能力。

#include "SF_systemValidator.h"
#include "solver/equation/SF_assemblyPlan.h"

#include <map>
#include <set>
#include <stdexcept>

namespace SF::System {

void validate(const ResolvedSimulationSystem& system) {
    std::set<std::string> unknowns;
    for (const auto& unknown : system.unknowns) {
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
    std::vector<std::string> equationIds;
    for (const auto& equation : system.equations) {
        if (equation.id.empty() || equation.name.empty()
            || !equations.insert(equation.id).second) {
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
    for (const auto& constraint : system.constraints) {
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
    for (const auto& block : system.solveBlocks) {
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
    // Build every plan during initialization.  A resolved equation may not be
    // explain-only metadata: it needs a stable lowering view and at least one
    // solve block that owns its execution.
    const Equation::AssemblyPlanRegistry plans(
        system.equationDefinitions,equationIds);
    for (const std::string& equation : equationIds) {
        (void)plans.at(equation);
        if (equationBlockReferences[equation] == 0) {
            throw std::runtime_error(
                "Resolved equation '"+equation
                +"' is not owned by any solve block.");
        }
    }
    for (const auto& requirement : system.requirements) {
        if (requirement.required && !requirement.available) {
            throw std::runtime_error(
                "Resolved system requires unavailable capability '"
                +requirement.name+"': "+requirement.detail+".");
        }
    }
    std::set<std::string> workspaces;
    for (const auto& workspace : system.workspaceRequirements) {
        if (workspace.id.empty() || workspace.components <= 0
            || !workspaces.insert(workspace.id).second) {
            throw std::runtime_error(
                "Resolved system has an invalid/duplicate workspace requirement.");
        }
    }
}

} // namespace SF::System
