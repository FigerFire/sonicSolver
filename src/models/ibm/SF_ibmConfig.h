/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.08.04-----------*/

#pragma once

/// @file SF_ibmConfig.h
/// @brief IBM/ILW 运行期只读配置，不依赖 parser 全局状态。

#include "SF_configTypes.h"

#include <algorithm>
#include <cmath>

/*
        storage:
        ┌──────── ghost / halo storage ────────┐
        │ NG                                   │
        │    ┌──── physical NX×NY×NZ ─────┐   │
        │    │                            │   │
        │    └────────────────────────────┘   │
        └─────────────────────────────────────┘
*/



namespace SF::IBM {

/// @brief IBM 几何分类、ILW 前处理和 ghost 重建共享的配置快照。
struct IBMRuntimeConfig {
    bool enabled = false;
    FDM::IBMMethod method = FDM::IBMMethod::Ghost;
    int requiredGhostLayers = 3;
    bool ilwEnabled = false;
    int ilwAccuracyOrder = 0;
    double gamma = 1.4;
    double normalAngleDegrees = 69.51268488527785;
    int normalSearchMinLayers = 1;
    int normalSearchMaxLayers = 0;
    int normalSearchTargetCandidates = 512;
    int normalSearchKeepSamples = 256;
    FDM::IBMForcingConfig forcing;

    /// @brief 返回 ILW Taylor 展开阶数；ILW3/5/7/9 对应 2/4/6/8。
    int requestedTaylorOrder() const {
        return ilwAccuracyOrder > 0
            ? std::max(1, ilwAccuracyOrder - 1)
            : 4;
    }

    /// @brief 把允许法向夹角转换为法向点积阈值。
    double normalAlignmentThreshold() const {
        constexpr double pi = 3.141592653589793238462643383279502884;
        return std::cos(normalAngleDegrees * pi / 180.0);
    }
};

} // namespace SF::IBM
