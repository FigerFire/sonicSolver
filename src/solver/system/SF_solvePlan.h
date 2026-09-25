#pragma once

/// @file SF_solvePlan.h
/// @brief Executable equations 与 execution policies 到 structured plan 的 lowering。

#include "core/system/SF_planFragment.h"
#include "SF_resolvedSimulationSystem.h"

namespace SF::System {

class SolvePlanner {
public:
    static std::vector<OpId> requiredOperations(const CompiledSolvePlan& plan);
    /// @brief 编译 plan。
    ///
    /// Planner 只做 merge/order/解析：片段（coupling/model 贡献）决定顺序与
    /// 重复次数，stage -> OpId 由 executable operation authority 解析。
    static CompiledSolvePlan compile(
        const ExecutableEquationSystem& system,
        const std::vector<ExecutionPolicy>& policies,
        const FDM::TimeRecipe& timeRecipe,
        const std::vector<PlanFragment>& fragments = {});
};

} // namespace SF::System
