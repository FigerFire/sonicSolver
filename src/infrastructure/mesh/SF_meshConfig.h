#pragma once

#include "SF_configTypes.h"

/// @file SF_meshConfig.h
/// @brief 网格加载、校验和分解使用的只读运行配置。

#include <array>
#include <string>
#include <vector>

namespace SF {

struct MeshRuntimeConfig {
    FDM::InitialConditionConfig initialConditions;
    double idealGasGamma = 1.4;
    double idealGasConstant = 287.05;
    int requiredGhostLayers = 1;
    std::string convectionScheme = "WENO5";
    int ilwOrder = 0;
    bool parallelEnabled = false;
    int parallelProcessCount = 1;
    bool automaticPartition = false;
    std::array<int, 3> partitionSplit{1, 1, 1};
    double haloTolerance = 1.0e-10;
    std::vector<std::string> configuredSetNames;
    std::vector<std::string> physicalBoundaryNames;
};

} // namespace SF
