/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.07.07-----------*/

#pragma once

/// @file SF_closure.h
/// @brief RPI 脱离直径、脱离频率和成核点密度闭式关系。

#include "SF_phaseChange.h"

#include <string>

namespace SF {
namespace Physics {
namespace PhaseChange {
namespace RPI {

/// @brief RPI 气泡脱离闭式参数。
struct BubbleClosure {
    double departureDiameter = 0.0;     ///< 脱离直径, m。
    double departureFrequency = 0.0;    ///< 脱离频率, 1/s。
    double nucleationSiteDensity = 0.0; ///< 成核点密度, 1/m2。
    double wallMassFlux = 0.0;          ///< 蒸发面质量通量, kg/(m2 s)。
};

/// @brief 规整 RPI 子模型选项。
std::string normalizeChoice(std::string value);

/// @brief 校验 RPI 闭式关系输入。
void validateClosureConfig(const Multiphase::PhaseChangeOptions& config,
                           const std::string& context);

/// @brief 按显式闭式关系计算局部气泡参数。
BubbleClosure evaluateClosure(const ModelContext& ctx,
                              double wallTemperature,
                              double liquidTemperature);

} // namespace RPI
} // namespace PhaseChange
} // namespace Physics
} // namespace SF
