/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_surfaceTension.h
/// @brief Level Set连续表面力(CSF)接口。

#include "SF_field.h"
#include "SF_state.h"

namespace SF {
namespace Physics {
namespace Multiphase {

/// @brief 表面张力体力工具。
class SurfaceTension {
public:
    /// @brief 计算一个点上的CSF体力 -sigma*kappa*delta(phi)*normal。
    /// @param levelSet 已更新几何的Level Set场。
    /// @param i,j,k 存储索引。
    /// @param sigma 表面张力系数, 必须非负。
    /// @param epsilon delta正则化半厚度, 必须为正。
    /// @return 单位体积表面张力体力。
    static Vector3 csfForce(const LevelSetField& levelSet,
                            int i, int j, int k,
                            double sigma,
                            double epsilon);
};

} // namespace Multiphase
} // namespace Physics
} // namespace SF
