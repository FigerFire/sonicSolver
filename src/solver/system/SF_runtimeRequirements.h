#pragma once

/// @file SF_runtimeRequirements.h
/// @brief RUNTIME — 编译/运行边界的 capability 结果、workspace 与
///        provider/runtime-service 需求。不属于 solve-plan IR。

#include "core/system/SF_equationIR.h"
#include "SF_numericalSystem.h"
#include "core/system/SF_solveProgram.h"

#include <algorithm>
#include <string>
#include <initializer_list>
#include <string_view>
#include <vector>

namespace SF::System {

enum class RuntimeStatus { Runnable, Unsupported, Invalid };
enum class BindingStatus { Resolved, Unsupported, Invalid };

/// @brief 编译期确定的 OpId -> provider 绑定；runtime 不重新选择实现。
struct ResolvedOperationBinding {
    OpId operation;
    std::string provider;
    BindingStatus status = BindingStatus::Unsupported;
    std::string reason;
};

struct RuntimeReport {
    RuntimeStatus status = RuntimeStatus::Unsupported;
    std::vector<OpId> requiredOperations;
    std::vector<OpId> missingOperations;
    std::string reason;
};

/// Diagnostic/validation summary; it never selects native runtime providers.
struct ExecutionCapabilitySignature {
    bool pressureConstraint = false;
    bool constantDensity = false;
    bool conservativeState = false;
    bool momentumPredictor = false;
    bool pressureCorrection = false;
    int outerCorrectors = 0;
    int pressureCorrectors = 0;
    int nonOrthogonalCorrectors = 0;
};

struct ExecutionRequirement {
    std::string name;
    bool required = false;
    bool available = false;
    std::string detail;
};

struct WorkspaceRequirement {
    std::string id;
    int components = 1;
    VariableLocation location = VariableLocation::EulerianCell;
    OwnershipKind ownership = OwnershipKind::EulerianGlobalDof;
};

/// @brief Concrete numerical resource required to bind the compiled program.
///
/// This is a small composition key, not a polymorphic "universal provider".
/// Mathematical content remains authoritative in ExecutableEquationSystem and
/// execution order remains authoritative in CompiledSolvePlan.
struct ProviderRequirement {
    std::string id;
    std::string responsibility;
};

/// @brief Runtime facility required by one or more concrete providers.
struct RuntimeServiceRequirement {
    std::string id;
    std::string reason;
};

/// @brief RUNTIME — 编译期产出的全部运行时需求。
///
/// 这一个值是 environment / execution / MPI backend 唯一需要读取的 RUNTIME
/// authority。它不持有 Field、mesh、communicator 或 case 输出设置。
struct RuntimeRequirements {
    RuntimeReport report;
    std::vector<ResolvedOperationBinding> operationBindings;
    ExecutionCapabilitySignature capabilities;
    std::vector<ExecutionRequirement> requirements;
    std::vector<WorkspaceRequirement> workspaceRequirements;
    std::vector<ProviderRequirement> providerRequirements;
    std::vector<RuntimeServiceRequirement> runtimeServiceRequirements;
};

/// @brief 四类编译产物的 typed 组合视图；consumer 只依赖自己需要的那一类。
///
/// 这里刻意持有引用而不是副本：CompiledSolvePlan 的 identity 是数值语义的一部分
/// （`bindSolvePlan` 校验收到的 plan 就是它的 plan），复制会破坏该不变量。
/// 禁止在此加入 CaseConfig、mesh、Field、communicator、IBM runtime state、
/// 输出设置、raw system、transformation history 或打印字符串。
struct Program {
    const ExecutableEquationSystem& equations;
    const CompiledNumericalSystem& numerics;
    const CompiledSolvePlan& solve;
    const RuntimeRequirements& runtime;
    const StateRegistry& state=equations.state;
};

/// @brief 只查询冻结计划实际请求的操作；不选择算法或补入生命周期。
inline bool requiresOperation(const RuntimeRequirements& requirements,std::string_view id) {
    const auto& operations=requirements.report.requiredOperations;
    return std::find(operations.begin(),operations.end(),id)!=operations.end();
}

/// @brief 已选 PressureMomentum 实现的资源契约；不按字段名或单个 OpId 启用算法。
/// 返回现有冻结数值输入；核验 STATE backing、完整操作组与持久面通量归属。
/// 没有该方法调用时返回 nullptr，不将未消费的默认配置变成活动资源。
const PressureNumericalConfig* requiredPressureInputs(const Program&);
/// @brief 从实际冻结相装配契约核验布局、旧时间层/通量与操作需求；不按配置启用。
const PhaseNumericalConfig* requiredPhaseInputs(const Program&);

/// @brief RUNTIME residency queries；只读取 RuntimeRequirements。
bool requiresCapability(const RuntimeRequirements& requirements,
                        std::string_view name);
bool requiresProvider(const RuntimeRequirements& requirements,
                      std::string_view id);
bool requiresRuntimeService(const RuntimeRequirements& requirements,
                            std::string_view id);

inline std::vector<OpId> assignedOperations(
        const RuntimeRequirements& requirements,
        std::initializer_list<std::string_view> providers) {
    std::vector<OpId> result;
    for (const auto& binding : requirements.operationBindings) {
        if (binding.status != BindingStatus::Resolved) continue;
        for (std::string_view provider : providers) {
            if (binding.provider == provider) {
                result.push_back(binding.operation);
                break;
            }
        }
    }
    return result;
}

} // namespace SF::System
