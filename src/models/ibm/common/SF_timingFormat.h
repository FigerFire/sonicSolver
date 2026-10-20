#pragma once

/// @file SF_timingFormat.h
/// @brief IBM setup 计时和秒值格式化小工具。

#include <chrono>
#include <iomanip>
#include <sstream>
#include <string>

namespace SF::IBM::Common {

inline double elapsedSeconds(
        std::chrono::steady_clock::time_point start) {
    return std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
}

inline std::string formatSeconds(double seconds) {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(6) << seconds << "s";
    return stream.str();
}

} // namespace SF::IBM::Common
