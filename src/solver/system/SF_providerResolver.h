#pragma once

/// @file SF_providerResolver.h
/// @brief Compiled operation capability 与现有 numerical provider 的唯一绑定点。

#include "SF_runtimeRequirements.h"
#include "SF_stateRealization.h"
#include "core/interfaces/SF_buildCapabilities.h"
#include "core/config/SF_equationComposition.h"

namespace SF::System {

/// Missing mathematical realizations fail before lowering unsupported operator operands.
std::string validateSolutionProviderContract(const StateRegistry& state,
    const EquationCompositionConfig& composition);

std::string validateCompiledConservativeStage(const ExecutableEquationSystem&,
    const CompiledSolvePlan&);

ExecutionCapabilitySignature compileExecutionCapabilities(
    const ExecutableEquationSystem& equations,
    const std::vector<LegacyExecutionPolicy>& policies,
    const CompiledSolvePlan* plan = nullptr);

/// Freeze final leaf ownership during compilation; runtime never calls this.
std::vector<ResolvedOperationBinding> compileOperationBindings(
    const ExecutableEquationSystem& equations,
    const CompiledNumericalSystem& numerics,
    CompiledSolvePlan& plan,
    const std::vector<LegacyExecutionPolicy>& policies,
    bool additionalContributions);

RuntimeReport reportOperationBindings(
    const CompiledSolvePlan& plan,
    const std::vector<ResolvedOperationBinding>& bindings);

} // namespace SF::System
