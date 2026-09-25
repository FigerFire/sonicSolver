/// @file SF_advect.h
/// @brief Level Set 输运、重初始化或界面几何模型实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

#include "SF_field.h"
#include "SF_state.h"

#include <vector>

namespace SF {
namespace Physics {
namespace Multiphase {

/// @brief Level Set 对流推进设置。
struct AdvectOptions {
    int order = 0;                ///< 用户输入的 HJ-WENO 阶数, 仅支持3/5/7。
    double wenoEpsilon = 0.0;     ///< 用户输入的 WENO 非线性权重正则量。
    double wenoPower = 0.0;       ///< 用户输入的 WENO 非线性权重指数。
    bool skipSolidCells = false;  ///< 用户输入：是否跳过 IBM solid/ghost 非流体点。
    bool skipSolidCellsDeclared = false; ///< 非流体点策略是否由外部配置传入。
};

/// @brief 装配 Hamilton–Jacobi WENO 空间 RHS。
///
/// 本类不拥有时间推进器。主流守恒量与 phi 必须由 solver 的同一显式 tableau
/// 更新，避免界面先走完整步、流场随后再走完整步的算子分裂。
class Advect {
public:
    /// @brief 用当前主流 stage 速度装配 d(phi)/dt。
    static void assembleRHS(const Field& field,
                            const LevelSetField& levelSet,
                            const AdvectOptions& options,
                            std::vector<double>& rhs);
};

} // namespace Multiphase
} // namespace Physics
} // namespace SF
