#pragma once

/// @file SF_workflow.h
/// @brief 求解 workflow 选择与能力兼容性检查。

#include "SF_config.h"
#include "SF_resolvedSimulationSystem.h"

#include <string>
#include <vector>

namespace SF::Workflow {

/// @brief 应用层可选的物理—方程 workflow（legacy，逐步被 stages 取代）。
/// @brief 一个已解析执行阶段的粗粒度类别。
enum class StageKind {
    Predictor,   ///< 显式守恒 predictor / RK 阶段
    Correction,  ///< pressure-velocity / PIMPLE 耦合校正
    Constraint,  ///< immersed / 代数约束求解
    Commit       ///< 最终状态提交
};

/// @brief 执行计划中的一个有序阶段：哪些方程/约束一起求解。
struct Stage {
    std::string id;
    StageKind kind = StageKind::Commit;
    std::vector<std::string> equations;
    std::vector<std::string> constraints;
};

/// @brief parser 和运行环境提供的中立 workflow 请求。
/// @brief 经能力检查后可供应用层执行的不可变计划。
struct Plan {
    /// @brief legacy dispatch key；runner 迁移到 stages 后移除。
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
/// @param request 中立物理/运行请求（legacy：仅用于推导 kind，runner 迁移后移除）。
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
