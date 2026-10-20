#pragma once

/// @file SF_solvePlan.h
/// @brief 已编译 method/temporal 控制流到唯一 structured plan。

#include "SF_resolvedSimulationSystem.h"
#include "SF_compiledTimeRecipe.h"

namespace SF::System {

class SolvePlanner {
public:
    static std::vector<OpId> requiredOperations(const CompiledSolvePlan& plan);
    /// @brief 编译 plan。
    ///
    /// HOW/WHICH 已决定顺序、重复次数和 OpId；planner 不推断 target 或旧策略。
    static CompiledSolvePlan compile(
        const ExecutableEquationSystem& system,
        const CompiledTimeRecipe& timeRecipe,
        const CompiledExecutionProgram& compiledExecution = {});
};

} // namespace SF::System
