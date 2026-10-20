/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.07.07-----------*/

#pragma once

/// @file SF_RPI.h
/// @brief RPI-like 壁面沸腾模型调度入口。

#include "SF_phaseChange.h"

#include <string>
#include <vector>

namespace SF {
namespace Physics {
namespace PhaseChange {
namespace RPI {

/// @brief 校验 RPI 模型配置。
/// @param config 相变配置。
/// @param context 报错上下文。
void validateConfig(const Multiphase::PhaseChangeOptions& config,
                    const std::string& context);

/// @brief 计算 RPI 壁面蒸发质量源。
/// @param ctx 相变模型上下文。
/// @param mdot 液相到气相质量源，单位 kg/(m3 s)。
void computeRates(const ModelContext& ctx, std::vector<double>& mdot);

} // namespace RPI
} // namespace PhaseChange
} // namespace Physics
} // namespace SF
