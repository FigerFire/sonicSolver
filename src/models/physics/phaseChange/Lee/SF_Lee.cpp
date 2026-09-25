/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.07.07-----------*/

/// @file SF_Lee.cpp
/// @brief Lee 相变模型调度实现。
///
/// addLeeRates 是外层的唯一调度入口。
/// 内部依次调用：
/// 1. SF_massTransfer::addMassSource   — 填充 mdot
///
/// 相质量更新只允许由 FluidStateModel 守恒源项统一提交。

#include "SF_Lee.h"

namespace SF {
namespace Physics {
namespace PhaseChange {

void addLeeRates(const ModelContext& ctx, std::vector<double>& mdot) {
    Lee::addMassSource(ctx, mdot);
}

} // namespace PhaseChange
} // namespace Physics
} // namespace SF
