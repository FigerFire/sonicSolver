/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_dimension.h
/// @brief 求解维度状态工具。

#include <array>
#include <cmath>
#include <string>

namespace SF {
namespace Math {

namespace Detail {
inline std::array<bool, 3>& activeDirections() {
    static std::array<bool, 3> state{true, true, true};
    return state;
}

inline std::array<double, 3>& inactiveDirectionNormal() {
    static std::array<double, 3> state{0.0, 0.0, 1.0};
    return state;
}
} // namespace Detail

/// @brief 重置为三维求解。
inline void resetActiveDirections() {
    Detail::activeDirections() = {true, true, true};
    Detail::inactiveDirectionNormal() = {0.0, 0.0, 1.0};
}

/// @brief 关闭一个坐标方向的空间离散。
inline void deactivateDirection(int axis) {
    if (axis >= 0 && axis < 3) {
        Detail::activeDirections()[(size_t)axis] = false;
    }
}

/// @brief 设置empty关闭方向在物理空间中的单位法向。
inline void setInactiveDirectionNormal(std::array<double, 3> normal) {
    const double mag = std::sqrt(normal[0] * normal[0]
                               + normal[1] * normal[1]
                               + normal[2] * normal[2]);
    if (mag <= 1.0e-14) {
        Detail::inactiveDirectionNormal() = {0.0, 0.0, 1.0};
        return;
    }
    Detail::inactiveDirectionNormal() = {
        normal[0] / mag, normal[1] / mag, normal[2] / mag};
}

/// @brief 返回empty关闭方向的物理空间单位法向。
inline std::array<double, 3> inactiveDirectionNormalVector() {
    return Detail::inactiveDirectionNormal();
}

/// @brief 查询坐标方向是否参与空间离散。
inline bool isDirectionActiveIndex(int axis) {
    return axis >= 0 && axis < 3
        && Detail::activeDirections()[(size_t)axis];
}

/// @brief 返回当前求解维度。
inline int activeDimensionCount() {
    int n = 0;
    for (bool active : Detail::activeDirections()) {
        if (active) ++n;
    }
    return n;
}

/// @brief 当且仅当当前为单一关闭方向的二维/准二维求解时，返回关闭方向。
/// @return 0=XI, 1=ETA, 2=ZETA；若不是恰好一个方向关闭则返回-1。
inline int singleInactiveDirectionIndex() {
    int inactive = -1;
    int count = 0;
    for (int axis = 0; axis < 3; ++axis) {
        if (!Detail::activeDirections()[(size_t)axis]) {
            inactive = axis;
            ++count;
        }
    }
    return count == 1 ? inactive : -1;
}

/// @brief 格式化当前active方向，用于日志。
inline std::string activeDirectionText() {
    std::string text;
    const auto& active = Detail::activeDirections();
    if (active[0]) text += text.empty() ? "XI" : "+XI";
    if (active[1]) text += text.empty() ? "ETA" : "+ETA";
    if (active[2]) text += text.empty() ? "ZETA" : "+ZETA";
    return text.empty() ? "None" : text;
}

} // namespace Math
} // namespace SF
