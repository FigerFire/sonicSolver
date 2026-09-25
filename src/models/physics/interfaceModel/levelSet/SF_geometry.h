/// @file SF_geometry.h
/// @brief Level Set 输运、重初始化或界面几何模型实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

#include "SF_field.h"
#include "SF_state.h"

namespace SF {
namespace Physics {
namespace Multiphase {

/// @brief Level Set 界面几何计算选项。
struct GeometryOptions {
    double gradientTolerance = 0.0; ///< 用户输入：界面邻域允许的最小 |grad(phi)|。
};

/// @brief 从 Level Set 场计算法向和曲率。
class Geometry {
public:
    /// @brief 计算 normal = grad(phi)/|grad(phi)| 和符号变化邻域的 kappa = div(normal)。
    /// @param field 参考网格。
    /// @param levelSet Level Set 场及几何缓存。
    /// @param options 几何计算选项。
    static void compute(const Field& field,
                        LevelSetField& levelSet,
                        const GeometryOptions& options);
};

} // namespace Multiphase
} // namespace Physics
} // namespace SF
