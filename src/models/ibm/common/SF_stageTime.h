#pragma once

/// @file SF_stageTime.h
/// @brief IBM stage 时间参数的显式合法性契约。

#include <cmath>
#include <stdexcept>

namespace SF::IBM::Common {

inline void validateStageTime(double time, double dt) {
    if (!std::isfinite(time) || !std::isfinite(dt) || dt < 0.0) {
        throw std::invalid_argument(
            "IBM stage requires finite time and non-negative dt.");
    }
}

} // namespace SF::IBM::Common
