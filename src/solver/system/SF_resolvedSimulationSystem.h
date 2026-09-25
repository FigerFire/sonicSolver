#pragma once

/// @file SF_resolvedSimulationSystem.h
/// @brief 启动阶段完整解析结果：每一阶段只有一个明确 authority。
///
/// WHAT / ORDER / HOW / RUNTIME 的分层 value objects 分别位于：
///   - SF_resolvedEquationSystem.h       (WHAT:  unknowns / equations / constraints)
///   - SF_solveProgram.h         (ORDER: policies / blocks / control-flow IR)
///   - SF_numericalSystem.h      (HOW:   term -> recipe bindings)
///   - SF_runtimeRequirements.h  (RUNTIME: capabilities / requirements)
///   - SF_transformation.h       (raw -> executable transformation IR)

#include "SF_caseClassification.h"
#include "SF_configTypes.h"
#include "SF_couplingStatus.h"
#include "core/system/SF_equationIR.h"
#include "SF_numericalSystem.h"
#include "SF_physicsTemplate.h"
#include "SF_runtimeRequirements.h"
#include "core/system/SF_solveProgram.h"
#include "SF_stateRealization.h"
#include "SF_transformation.h"

#include <string>
#include <string_view>
#include <vector>

namespace SF::System {

struct ResolvedSimulationSystem {
    FDM::TimeRecipe timeRecipe;
    /// @brief STATE REALIZATION：由 executable equations 的 role 导出。
    ///        它不是 solver family，也不选择运行时路径。
    CompiledStateRealization realization;
    /// @brief 压力耦合 preset 注册状态（active/inactive/invalid/unsupported）。
    CouplingReport coupling;
    /// EXPLAIN-ONLY：不参与执行，也不在 Program 中可见。
    CaseClassification classification;
    RawEquationSystem rawSystem;
    ExecutableEquationSystem executableSystem;
    CompiledNumericalSystem numericalSystem;
    std::vector<TransformationRecord> transformations;
    std::vector<ExecutionPolicy> executionPolicies;
    CompiledSolvePlan solvePlan;
    RuntimeRequirements runtime;
};

/// @brief §19 编译产物视图：raw system + transformation history + typed program。
///
/// 与 Program 一样持有引用，只作为 builder 输出的分层观察入口，
/// 不复制任何编译结果，也不重新引入第二份 authority。
struct CompilationResult {
    const RawEquationSystem& raw;
    const std::vector<TransformationRecord>& transformations;
    Program program;
};

inline Program programOf(const ResolvedSimulationSystem& system) {
    return Program{system.executableSystem, system.numericalSystem,
                   system.solvePlan, system.runtime};
}

inline CompilationResult compilationOf(const ResolvedSimulationSystem& system) {
    return CompilationResult{system.rawSystem, system.transformations,
                             programOf(system)};
}

bool hasUnknown(const ResolvedSimulationSystem& system, std::string_view id);
bool hasEquation(const ResolvedSimulationSystem& system, std::string_view id);
const Equation::Definition& equationDefinition(
    const ResolvedSimulationSystem& system, std::string_view id);
bool hasEquationPrefix(
    const ResolvedSimulationSystem& system, std::string_view prefix);
bool hasConstraint(const ResolvedSimulationSystem& system, std::string_view id);
bool requiresCapability(
    const ResolvedSimulationSystem& system, std::string_view name);
bool requiresProvider(
    const ResolvedSimulationSystem& system, std::string_view id);
bool requiresRuntimeService(
    const ResolvedSimulationSystem& system, std::string_view id);

} // namespace SF::System
