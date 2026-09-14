#pragma once

/// @file SF_workflow.h
/// @brief 求解 workflow 选择与能力兼容性检查。

#include "SF_config.h"

#include <string>

namespace SF::Workflow {

/// @brief 应用层可选的物理—方程 workflow。
enum class Kind {
    SingleFluid,
    OneFluidInterface,
    LegacyMultiphase,
    Homogeneous,
    EulerianEulerian
};

/// @brief parser 和运行环境提供的中立 workflow 请求。
struct Request {
    bool multiphaseEnabled = false;
    bool homogeneous = false;
    bool eulerianEulerian = false;
    bool interfaceModel = false;
    bool interfaceGhostFluid = false;
    bool legacyMixture = false;
    bool phaseChange = false;
    bool transportedLegacyAlpha = false;
    bool ibm = false;
    bool ibmForcing = false;
    bool ilw = false;
};

/// @brief 经能力检查后可供应用层执行的不可变计划。
struct Plan {
    Kind kind = Kind::SingleFluid;
    bool transportedLegacyAlpha = false;
    bool turbulenceHasTransportState = false;

    /// @brief 返回稳定的日志名称。
    std::string name() const;
};

/// @brief 在网格和求解状态分配前选择并验证 workflow。
/// @param request 中立物理/运行请求。
/// @param config 已物化的求解器配置。
/// @return 通过基本能力检查的 workflow 计划。
/// @throws std::runtime_error 不支持的组合不会被替换为其他 workflow。
Plan makePlan(const Request& request, const FDM::SolverConfig& config);

/// @brief 用实际网格执行形态补充验证 workflow。
/// @param plan 已选择的 workflow。
/// @param multiField 同一执行分支是否持有多个结构 patch。
/// @param hasCoupledInterfaces 网格是否含跨 patch 耦合接口。
/// @throws std::runtime_error 当前 workflow 不支持该执行形态时抛出。
void validateMeshExecution(
    const Plan& plan,
    bool multiField,
    bool hasCoupledInterfaces);

} // namespace SF::Workflow
