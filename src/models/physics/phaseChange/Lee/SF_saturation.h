/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.07.07-----------*/

/// @file SF_saturation.h
/// @brief Lee 模型的饱和物性计算模块。
///
/// 使用 Clausius-Clapeyron 关系从局部压力计算饱和温度。
/// 若未配置饱和压力，则直接返回配置中的常饱和温度。

#pragma once

#include "SF_phaseChange.h"

#include <cmath>
#include <stdexcept>

namespace SF {
namespace Physics {
namespace PhaseChange {
namespace Lee {

/// @brief Clausius-Clapeyron 饱和温度计算。
///
/// 从局部压力 p 和参考饱和点 (psat0, Tsat0) 计算 Tsat(p)：
///   Tsat = 1 / (1/Tsat0 - R * ln(p / psat0) / L)
///
/// @param p 局部压力。
/// @param psat0 参考饱和压力。
/// @param tsat0 参考饱和温度 (K)。
/// @param L 潜热 (J/kg)。
/// @param R 气体常数 (J/(kg*K)), 默认 461.52 (水蒸气)。
/// @return 局部饱和温度 (K)。
inline double tsatFromClausiusClapeyron(double p, double psat0,
                                        double tsat0, double L,
                                        double R = 461.52) {
    if (p <= 0.0 || psat0 <= 0.0 || !std::isfinite(p)) {
        return tsat0;
    }
    if (!std::isfinite(psat0) || !std::isfinite(tsat0) ||
        !std::isfinite(L) || !std::isfinite(R)) {
        throw std::runtime_error(
            "Lee::tsatFromClausiusClapeyron: non-finite input.");
    }
    if (tsat0 <= 0.0 || L <= 0.0 || R <= 0.0) {
        throw std::runtime_error(
            "Lee::tsatFromClausiusClapeyron: Tsat0, L, R must be > 0.");
    }
    const double ratio = p / psat0;
    if (ratio <= 0.0) return tsat0;
    const double logTerm = std::log(ratio) * R / L;
    const double invTsat = 1.0 / tsat0 - logTerm;
    if (invTsat <= 0.0) {
        return tsat0;
    }
    const double tsat = 1.0 / invTsat;
    if (!std::isfinite(tsat) || tsat <= 0.0) return tsat0;
    return tsat;
}

/// @brief 根据配置获取饱和温度。
///
/// 当 config 提供了 saturationPressure > 0 时启用 Clausius-Clapeyron 修正，
/// 否则直接返回常饱和温度。
///
/// @param p 局部压力。
/// @param pc PhaseChangeOptions 配置。
/// @param gasConstant 可选的显式气体常数；0 时使用默认值 461.52。
/// @return 饱和温度 (K)。
inline double saturationTemperature(double p,
                                    const Multiphase::PhaseChangeOptions& pc,
                                    double gasConstant = 0.0) {
    if (!std::isfinite(pc.saturationTemperature) ||
        pc.saturationTemperature <= 0.0) {
        throw std::runtime_error(
            "Lee::saturationTemperature: invalid Tsat in config.");
    }
    if (pc.saturationPressure > 0.0 && pc.latentHeat > 0.0) {
        const double R = (gasConstant > 0.0) ? gasConstant : 461.52;
        return tsatFromClausiusClapeyron(p, pc.saturationPressure,
                                         pc.saturationTemperature,
                                         pc.latentHeat, R);
    }
    return pc.saturationTemperature;
}

} // namespace Lee
} // namespace PhaseChange
} // namespace Physics
} // namespace SF
