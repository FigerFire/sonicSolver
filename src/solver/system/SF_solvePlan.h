#pragma once

/// @file SF_solvePlan.h
/// @brief Executable equations 与 execution policies 到 structured plan 的 lowering。

#include "SF_resolvedSimulationSystem.h"

namespace SF::System {

class SolvePlanner {
public:
    static std::vector<OpId> requiredOperations(const CompiledSolvePlan& plan);
    static ExecutionCapabilitySignature capabilities(
        const ExecutableEquationSystem& system,
        const std::vector<ExecutionPolicy>& policies);

    static CompiledSolvePlan compile(
        const ExecutableEquationSystem& system,
        const std::vector<ExecutionPolicy>& policies,
        std::string_view timeIntegrator);
};

} // namespace SF::System
