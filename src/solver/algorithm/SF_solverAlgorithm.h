#pragma once

/// @file SF_solverAlgorithm.h
/// @brief 流动求解族算法工厂调度入口。

#include "SF_config.h"
#include "SF_flowAlgorithm.h"
#include "solver/algorithm/densityBased/SF_densityBased.h"
#include "solver/algorithm/pressureBased/SF_pressureBased.h"
#include "solver/algorithm/SF_compressible.h"

#include <memory>

namespace SF::FDM {

/// @brief 按强类型配置创建密度基或压力基算法对象。
std::unique_ptr<IFlowAlgorithm> makeFlowAlgorithm(
    const SolverConfig& config);

} // namespace SF::FDM
