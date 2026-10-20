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
    CompiledSolvePlan result;
    SolveBlock time;
    time.id="TIME_RECIPE";
    time.name=std::string(FDM::toString(timeRecipe.id()))+" time recipe";
    time.strategy=FDM::toString(timeRecipe.id());
    time.strategyKind=FDM::SolveStrategyKind::ExplicitTimeIntegration;
    for (const auto& step:compiledExecution.steps) {
        for (const auto& input:step.sourceMathInputs)
            if (system.registry.contains(input)) time.equations.push_back(input);
        for (const auto& call:step.calls)
            if (system.registry.contains(call.equation)) time.equations.push_back(call.equation);
        if (system.state.contains(step.source.target.symbol))
            time.unknowns.push_back(step.source.target.symbol);
    }
    result.blocks.push_back(std::move(time));
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
