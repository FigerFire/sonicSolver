#pragma once

/// @file SF_runtimeConfig.h
/// @brief 将 CaseConfig 与 solver configuration 转换为 mesh/IBM runtime 配置。
///
/// Data flow:
///   CaseConfig + SolverConfig
///       -> MeshRuntimeConfig / IBMRuntimeConfig
///       -> application environment construction
///
/// 本文件不执行 solver、不选择 runner，也不生成日志或输出字段。

#include "SF_config.h"
#include "SF_meshConfig.h"
#include "SF_ibmConfig.h"
#include "SF_caseConfig.h"

namespace SF::Application::Runtime {

MeshRuntimeConfig makeMeshRuntimeConfig(
    const CaseConfig& caseConfig,
    const FDM::SolverConfig& solverConfig,
    int requiredTermHaloWidth);
IBM::IBMRuntimeConfig makeIBMRuntimeConfig(
    const FDM::SolverConfig& solverConfig,
    int requiredTermHaloWidth);

} // namespace SF::Application::Runtime
