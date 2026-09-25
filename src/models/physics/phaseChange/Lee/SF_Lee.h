/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.07.07-----------*/

/// @file SF_Lee.h
/// @brief Lee 相变模型调度入口。
///
/// 本文件是 Lee 文件夹的对外调度接口，聚合以下子模块：
/// - SF_saturation.h   : 饱和物性 (Clausius-Clapeyron)
/// - SF_massTransfer.h : 质量传递速率
/// 潜热由 EOS referenceEnergy 表达，不存在独立 thermal source。
///
/// 所有离散和算子调用来自 numerics/ 和 math/ 模块。

#pragma once

#include "SF_phaseChange.h"
#include "SF_saturation.h"
#include "SF_massTransfer.h"

#include <vector>

namespace SF {
namespace Physics {
namespace PhaseChange {

/// @brief Lee 模型主入口：填充 mdot 质量源数组。
///
/// @param ctx 相变模型上下文。
/// @param mdot 输出的质量源项数组 (kg/(m^3*s))。
void addLeeRates(const ModelContext& ctx, std::vector<double>& mdot);

} // namespace PhaseChange
} // namespace Physics
} // namespace SF
