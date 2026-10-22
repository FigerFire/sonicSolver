/// @file SF_solvePlan.cpp
/// @brief 冻结已编译 method/temporal 控制流，不推断旧策略或数学 target。
#include "SF_solvePlan.h"
#include "SF_methodObjects.h"
#include <algorithm>
#include <stdexcept>
namespace SF::System {
CompiledSolvePlan SolvePlanner::compile(
        const ExecutableEquationSystem& system,
        const CompiledTimeRecipe& timeRecipe,
        const CompiledExecutionProgram& compiledExecution) {
    if (compiledExecution.steps.empty())
        throw std::runtime_error("SolvePlanner requires compiled HOW equation occurrences; no legacy lowering is available.");
    if (timeRecipe.stageCount()<1) throw std::runtime_error("Invalid configured time recipe.");
    CompiledSolvePlan result;
    result.root=compiledExecution.hasTemporalRoot
        ? compiledExecution.temporalRoot : compileMethodProgram(compiledExecution);
    return result;
}

std::vector<OpId> SolvePlanner::requiredOperations(const CompiledSolvePlan& plan) {
    std::vector<OpId> result;
    const auto visit = [&](const auto& self, const SolvePlanNode& node) -> void {
        if (!node.operation.empty()
            && std::find(result.begin(),result.end(),node.operation) == result.end()) {
            result.push_back(node.operation);
        }
        for (const auto& child : node.children) self(self,child);
    };
    visit(visit,plan.root);
    return result;
}

} // namespace SF::System
