#pragma once

/// @file SF_providerResolver.h
/// @brief Compiled operation capability 与现有 numerical provider 的唯一绑定点。

#include "SF_runtimeRequirements.h"
#include "SF_stateRealization.h"

namespace SF::System {

ExecutionCapabilitySignature compileExecutionCapabilities(
    const ExecutableEquationSystem& equations,
    const std::vector<ExecutionPolicy>& policies);

std::vector<ResolvedOperationBinding> resolveOperationBindings(
    const ExecutableEquationSystem& equations,
    const CompiledStateRealization& realization,
    const CompiledNumericalSystem& numerics,
    const CompiledSolvePlan& plan,
    const std::vector<ExecutionPolicy>& policies);

RuntimeReport reportOperationBindings(
    const CompiledSolvePlan& plan,
    const std::vector<ResolvedOperationBinding>& bindings);

} // namespace SF::System
