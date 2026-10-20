#pragma once

/// @file SF_report.h
/// @brief 把已解析的 solver/runtime 配置格式化为人类可读运行信息。
///
/// Data flow:
///   CaseConfig / SolverConfig / timestep values
///       -> stable text formatting
///       -> CLI and runtime log
///
/// 本文件只观察配置和值，不修改 state、控制 timestep 或选择执行路径。

#include "SF_config.h"
#include "SF_caseConfig.h"

#include <string>

namespace SF {
namespace Physics::Multiphase { struct MultiPhaseConfig; }
namespace Application::Report {

/// @brief 多相模型摘要；flowLabel 是兼容输入标签（可为空），不是求解器身份。
std::string multiPhaseSummary(
    const Physics::Multiphase::MultiPhaseConfig& config,
    const std::string& flowLabel);
std::string formatTimeValue(double value);
std::string formatTimeStepStatus(double time, double dt);
void broadcastSolverConfig(
    const FDM::SolverConfig& config,
    bool eulerianEulerian = false);
std::string formatRunControl(const CaseConfig& config);

} // namespace Application::Report
} // namespace SF
