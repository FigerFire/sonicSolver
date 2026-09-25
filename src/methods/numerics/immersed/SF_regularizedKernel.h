#pragma once

/// @file SF_regularizedKernel.h
/// @brief IBM 紧支撑核与体积加权归一化的无状态数值算子。

#include <cmath>
#include <stdexcept>

namespace SF::FDM::Immersed {

/// @brief Wendland C2 紧支撑径向核；输入和支撑半径使用相同长度单位。
inline double wendlandC2(double distance, double supportRadius) {
    if (!std::isfinite(distance) || distance < 0.0
        || !std::isfinite(supportRadius) || supportRadius <= 0.0) {
        throw std::runtime_error(
            "IBM Wendland C2 kernel requires finite distance and positive "
            "support radius.");
    }
    const double q = distance/supportRadius;
    if (q >= 1.0) return 0.0;
    const double oneMinus = 1.0-q;
    return oneMinus*oneMinus*oneMinus*oneMinus*(1.0+4.0*q);
}

/// @brief 构造插值行的未归一化 `delta_h(x-X) V_i` 权重。
inline double volumeWeightedKernel(double kernelValue, double dualVolume) {
    if (!std::isfinite(kernelValue) || kernelValue < 0.0
        || !std::isfinite(dualVolume) || dualVolume <= 0.0) {
        throw std::runtime_error(
            "IBM kernel weighting requires non-negative kernel and positive "
            "Eulerian dual volume.");
    }
    return kernelValue*dualVolume;
}

} // namespace SF::FDM::Immersed
