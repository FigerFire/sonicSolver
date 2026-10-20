/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_weightTypes.h
/// @brief IBM ghost-cell 权重和分类数据结构。

#include "geoProcessing/SF_STLGeometry.h"

#include <vector>

namespace SF {
namespace IBM {
namespace GhostIBM {

using GeoProcessing::Point;

/// @brief 一个流体 donor 单元及其插值权重。
struct DonorWeight {
    int i = 0;
    int j = 0;
    int k = 0;
    double weight = 0.0;
};

/// @brief 一个 IBM ghost-cell 的几何截距、镜像点和 donor 权重。
struct GhostCellWeight {
    int ghostI = 0;
    int ghostJ = 0;
    int ghostK = 0;
    Point imagePoint;
    Point boundaryIntercept;
    Point wallNormal;
    double wallDistance = 0.0;
    std::vector<DonorWeight> donors;
};

/// @brief IBM几何分类计数。
struct ClassificationCounts {
    int fluid = 0;
    int ghost = 0;
    int solid = 0;
};

} // namespace GhostIBM
} // namespace IBM
} // namespace SF
