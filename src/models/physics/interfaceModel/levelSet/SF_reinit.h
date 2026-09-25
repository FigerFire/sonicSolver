/// @file SF_reinit.h
/// @brief Level Set 输运、重初始化或界面几何模型实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

#include "SF_field.h"
#include "SF_state.h"

#include <functional>

namespace SF {
namespace Physics {
namespace Multiphase {

/// @brief Level Set 重初始化设置。
struct ReinitOptions {
    int pseudoSteps = 0;          ///< 伪时间步数, 必须大于0。
    double pseudoTimeStep = 0.0;  ///< 伪时间步长, 必须为正。
    int order = 0;                ///< 用户输入的 HJ-WENO 阶数, 仅支持3/5/7。
    double wenoEpsilon = 0.0;     ///< 用户输入的 WENO 非线性权重正则量。
    double wenoPower = 0.0;       ///< 用户输入的 WENO 非线性权重指数。
    double signSmoothingFactor = 0.0; ///< S(phi0) 平滑宽度相对最小网格尺度的倍率。
    std::function<void()> prepareStage; ///< 每个伪时间步的边界与halo同步。
};

/// @brief Sussman/Fatemi 风格 signed-distance 重初始化。
class Reinit {
public:
    /// @brief 求解 phi_tau + S(phi0)(|grad(phi)| - 1) = 0。
    /// @param field 参考网格。
    /// @param levelSet Level Set 场。
    /// @param options 重初始化设置。
    static void advance(const Field& field,
                        LevelSetField& levelSet,
                        const ReinitOptions& options);
};

} // namespace Multiphase
} // namespace Physics
} // namespace SF
