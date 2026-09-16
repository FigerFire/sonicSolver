#pragma once

/// @file SF_workflow.h
/// @brief 求解 workflow 选择与能力兼容性检查。

#include "SF_config.h"
#include "SF_resolvedSimulationSystem.h"
#include "core/interfaces/SF_interfaces.h"

#include <string>
#include <vector>

namespace SF::Workflow {

/// @brief Workflow 与 runtime stepper 共用的 solve-stage 表示。
using StageKind = FDM::SolveStageKind;
using Stage = FDM::SolveStage;

/// @brief parser 和运行环境提供的中立 workflow 请求。
/// @brief 经能力检查后可供应用层执行的不可变计划。
struct Plan {
    /// @brief 由 ResolvedSimulationSystem 解析出的有序执行阶段。
    std::vector<Stage> stages;
    /// @brief 来自 ResolvedSimulationSystem 的时间积分器名称。
    std::string timeIntegrator;

    /// @brief 返回稳定的日志名称。
    std::string name() const;
};

/// @brief 由已解析的数学系统生成执行计划，并完成能力检查。
/// @param system 已冻结的 ResolvedSimulationSystem。
/// @param config 已物化的求解器配置。
/// @return 通过基本能力检查的 workflow 计划。
/// @throws std::runtime_error 不支持的组合不会被替换为其他 workflow。
Plan makePlan(
    const System::ResolvedSimulationSystem& system,
    const FDM::SolverConfig& config);

/// @brief 用实际网格执行形态补充验证 workflow。
/// @param plan 已选择的 workflow。
/// @param multiField 同一执行分支是否持有多个结构 patch。
/// @param hasCoupledInterfaces 网格是否含跨 patch 耦合接口。
/// @throws std::runtime_error 当前 workflow 不支持该执行形态时抛出。
void validateMeshExecution(
    const System::ResolvedSimulationSystem& system,
    const FDM::SolverConfig& config,
    bool multiField,
    bool hasCoupledInterfaces);

} // namespace SF::Workflow
