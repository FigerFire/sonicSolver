/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.07.07-----------*/

#pragma once

/// @file SF_wallMapping.h
/// @brief RPI 壁面 patch 到近壁流体点的拓扑映射。

#include "SF_phaseChange.h"

#include <vector>

namespace SF {
namespace Physics {
namespace PhaseChange {
namespace RPI {

/// @brief 一个壁面点及其相邻求解点。
struct WallCellSample {
    int wallI = 0, wallJ = 0, wallK = 0;
    int fluidI = 0, fluidJ = 0, fluidK = 0;
    double wallNormalDistance = 0.0;
    double areaOverVolume = 0.0; ///< 壁面面积/相邻控制体体积，单位 1/m。
};

/// @brief 将壁面 set 点映射到相邻流体点。
std::vector<WallCellSample> mapWallPatch(const ModelContext& ctx,
                                         const WallHeatSetting& setting);

} // namespace RPI
} // namespace PhaseChange
} // namespace Physics
} // namespace SF
