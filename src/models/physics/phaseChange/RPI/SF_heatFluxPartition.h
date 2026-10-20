/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.07.07-----------*/

#pragma once

/// @file SF_heatFluxPartition.h
/// @brief RPI 壁面热流分配。

#include "SF_closure.h"

#include <string>

namespace SF {
namespace Physics {
namespace PhaseChange {
namespace RPI {

/// @brief RPI 三部分热流。
struct HeatFluxPartition {
    double convective = 0.0;   ///< 单相对流热流 q_c, W/m2。
    double quenching = 0.0;    ///< 淬冷热流 q_q, W/m2。
    double evaporative = 0.0;  ///< 蒸发潜热热流 q_e, W/m2。
    double total = 0.0;        ///< q_c + q_q + q_e。
};

/// @brief 恒热流 RPI 壁面闭合后的温度、气泡参数和热流账本。
struct WallHeatBalance {
    double wallTemperature = 0.0; ///< 平衡后的局部壁温, K。
    BubbleClosure closure;        ///< 平衡壁温下的气泡闭式量。
    HeatFluxPartition heat;       ///< 平衡后的三部分热流。
};

/// @brief 校验热流分配模型输入。
void validateHeatFluxConfig(const Multiphase::PhaseChangeOptions& config,
                            const std::string& context);

/// @brief 计算 RPI 热流分配。
HeatFluxPartition evaluateHeatFluxPartition(const ModelContext& ctx,
                                            double wallTemperature,
                                            double liquidTemperature,
                                            double wallMassFlux);

/// @brief 按 wallTemperatureModel 计算显式壁温或恒热流平衡壁温。
///
/// heatFluxBalance 使用有界二分求解
/// q_conv(Tw)+q_quench(Tw)+q_evap(Tw)=q_wall；无物理解时 fail-fast。
WallHeatBalance evaluateWallHeatBalance(const ModelContext& ctx,
                                        double liquidTemperature,
                                        double referenceWallTemperature,
                                        double wallHeatFlux);

/// @brief 返回是否执行总热流预算检查。
bool budgetCheckEnabled(const Multiphase::PhaseChangeOptions& config);

} // namespace RPI
} // namespace PhaseChange
} // namespace Physics
} // namespace SF
