#pragma once

/// @file SF_numericsPolicy.h
/// @brief 仅由输入参数决定的数值模板与 ILW 阶数规则。

#include "core/config/SF_configTypes.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace SF::FDM {

/// @brief 判断 ILW 精度阶是否受支持。
inline bool isSupportedILWOrder(int order) {
    return order == 0 || order == 3 || order == 5
        || order == 7 || order == 9;
}

/// @brief 将 ILW 精度阶映射成 Taylor 导数最高阶。
inline int ilwTaylorOrderFromAccuracy(int order) {
    if (order <= 0) return 0;
    return std::max(1, order - 1);
}

/// @brief 根据 ILW 精度返回 Taylor ghost 外推需要的最小虚胞层数。
inline int requiredGhostLayersForILW(int order) {
    return order <= 0 ? 0 : (order + 1) / 2;
}

/// @brief Fail fast for declared numerical modes that do not have an implementation.
/// @param numerics Parsed numerical-method configuration.
inline void validateNumericsConfig(const NumericsConfig& numerics) {
    if (!std::isfinite(numerics.cfl) || numerics.cfl <= 0.0) {
        throw std::invalid_argument(
            "CFL/maxCo must be a finite positive value.");
    }
    if (!std::isfinite(numerics.maxDeltaT)
        || numerics.maxDeltaT <= 0.0) {
        throw std::invalid_argument(
            "maxDeltaT must be a finite positive value.");
    }
    if (!std::isfinite(numerics.idealGasGamma)
        || numerics.idealGasGamma <= 1.0
        || !std::isfinite(numerics.idealGasConstant)
        || numerics.idealGasConstant <= 0.0) {
        throw std::invalid_argument(
            "legacy ideal-gas numerics require gamma > 1 and R > 0.");
    }
    if (numerics.formulation
        != EquationFormulation::ConservativeFluxDifference) {
        throw std::invalid_argument(
            "primitive/non-conservative differential formulation is not "
            "supported for the density-based compressible solver; use "
            "conservativeFluxDifference.");
    }
    if (numerics.reconstruction
        != ReconstructionVariable::Characteristic) {
        throw std::invalid_argument(
            "requested face reconstruction '"
            + std::string(toString(numerics.reconstruction))
            + "' is not implemented yet; current WENO flux assembly supports "
            "characteristic reconstruction.");
    }
}

} // namespace SF::FDM
