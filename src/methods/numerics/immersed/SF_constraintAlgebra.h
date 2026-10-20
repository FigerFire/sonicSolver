#pragma once

/// @file SF_constraintAlgebra.h
/// @brief IBM 局部 Schur 补和约束质量代数的无状态算子。

#include <cmath>
#include <stdexcept>

namespace SF::FDM::Immersed {

/// @brief 计算局部近似 `J M^-1 J^T M_L` 的质量对角。
template <typename Row, typename InverseMassGetter>
double massDiagonal(const Row& row, double measure,
                    InverseMassGetter&& inverseMass) {
    if (!std::isfinite(measure) || measure<=0.0) {
        throw std::runtime_error(
            "IBM mass diagonal requires a positive marker measure.");
    }
    double diagonal=0.0;
    for (const auto& entry:row) {
        diagonal+=entry.value*entry.value*measure
                 *inverseMass(entry.cell);
    }
    if (!std::isfinite(diagonal) || diagonal<=0.0) {
        throw std::runtime_error(
            "IBM interpolation produced a singular mass diagonal.");
    }
    return diagonal;
}

} // namespace SF::FDM::Immersed
