#pragma once

/// @file SF_inspection.h
/// @brief application 层统一解析 case 的检查接口。

#include "SF_caseConfig.h"
#include "SF_ibmConfig.h"
#include "SF_ibmExplain.h"
#include "SF_meshConfig.h"
#include "SF_resolvedSimulationSystem.h"

#include <string>

namespace SF::Application {

/// @brief `check`、`explain`、`doctor` 和 `run` 共享的解析快照。
///
/// 该对象只包含启动阶段的 value objects；不会持有 Field、MPI communicator
/// 或 parser 全局状态，因此 CLI、GUI 和实际求解可以共享同一套解析结果。
struct CaseInspection {
    CaseConfig config;
    MeshRuntimeConfig mesh;
    IBM::IBMRuntimeConfig ibm;
    IBM::ImmersedExplain ibmExplain;
    System::ResolvedSimulationSystem system;
};

/// @brief 从已经读取的 CaseConfig 解析并校验完整执行系统。
/// @throws std::exception 输入、workflow 或能力组合不支持时抛出。
CaseInspection inspectCase(const CaseConfig& config);

/// @brief 加载并解析 case；native YAML 和 compatibility case 共用此入口。
/// @param path case 目录或兼容格式入口文件。
/// @throws std::exception case 无法读取或组合不支持时抛出。
CaseInspection inspectCase(const std::string& path);

} // namespace SF::Application
