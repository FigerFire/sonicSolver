#pragma once

/// @file SF_interpolation.h
/// @brief IBM 稀疏插值 `J u` 的无状态数值算子。

#include "SF_valueTypes.h"

#include <cmath>
#include <stdexcept>

namespace SF::FDM::Immersed {

/// @brief 计算 Lagrangian 点速度 `J u`。
template <typename Row, typename VelocityGetter>
Vector3 interpolate(const Row& row, VelocityGetter&& velocity) {
    Vector3 result;
    for (const auto& entry:row) {
        if (!std::isfinite(entry.value)) {
            throw std::runtime_error(
                "IBM interpolation row contains a non-finite weight.");
        }
        result=result+velocity(entry.cell)*entry.value;
    }
    return result;
}

} // namespace SF::FDM::Immersed
