#pragma once

/// @file SF_spreading.h
/// @brief IBM 质量加权伴随传播 `J^T M_L Lambda` 的无状态算子。

#include "SF_valueTypes.h"

#include <cmath>
#include <stdexcept>

namespace SF::FDM::Immersed {

/// @brief 产生积分力贡献；调用者决定串行累加或 Runtime sum-to-owner。
template <typename Row, typename Accumulator>
void spread(const Row& row, double measure, const Vector3& multiplier,
            Accumulator&& accumulate) {
    if (!std::isfinite(measure) || measure<=0.0) {
        throw std::runtime_error(
            "IBM spreading requires a positive marker measure.");
    }
    for (const auto& entry:row) {
        accumulate(entry.cell,multiplier*(measure*entry.value));
    }
}

} // namespace SF::FDM::Immersed
