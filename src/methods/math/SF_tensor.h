#pragma once

/// @file SF_tensor.h
/// @brief 三维张量的纯代数算子；不依赖 Field、Mesh 或通信模块。

#include <array>
#include <cstddef>

#include "frame/SF_vector.h"

namespace SF::Math {

/// @brief 三维二阶张量；第一下标为分量，第二下标为微分方向。
struct Tensor3 {
    std::array<Vector3, 3> rows{};

    Vector3& operator[](std::size_t row) { return rows[row]; }
    const Vector3& operator[](std::size_t row) const { return rows[row]; }
};

/// @brief 三维对称二阶张量，按六个独立分量存储。
struct SymmTensor3 {
    double xx = 0.0;
    double yy = 0.0;
    double zz = 0.0;
    double xy = 0.0;
    double xz = 0.0;
    double yz = 0.0;
};

/// @brief 返回二阶张量的转置。
inline Tensor3 transpose(const Tensor3& value) {
    Tensor3 result;
    for (std::size_t row = 0; row < 3; ++row) {
        for (std::size_t column = 0; column < 3; ++column) {
            result[row][column] = value[column][row];
        }
    }
    return result;
}

/// @brief 返回 `0.5*(A+A^T)`。
inline SymmTensor3 symm(const Tensor3& value) {
    return {
        value[0][0], value[1][1], value[2][2],
        0.5 * (value[0][1] + value[1][0]),
        0.5 * (value[0][2] + value[2][0]),
        0.5 * (value[1][2] + value[2][1])};
}

/// @brief 返回 `A+A^T`；名称显式保留二倍系数。
inline SymmTensor3 twoSymm(const Tensor3& value) {
    return {
        2.0 * value[0][0],
        2.0 * value[1][1],
        2.0 * value[2][2],
        value[0][1] + value[1][0],
        value[0][2] + value[2][0],
        value[1][2] + value[2][1]};
}

/// @brief 返回二阶张量迹。
inline double trace(const Tensor3& value) {
    return value[0][0] + value[1][1] + value[2][2];
}

/// @brief 返回对称张量迹。
inline double trace(const SymmTensor3& value) {
    return value.xx + value.yy + value.zz;
}

/// @brief 对称张量双点积 `A:A`，非对角项计两次。
inline double doubleDot(const SymmTensor3& first,
                        const SymmTensor3& second) {
    return first.xx * second.xx
         + first.yy * second.yy
         + first.zz * second.zz
         + 2.0 * (first.xy * second.xy
                + first.xz * second.xz
                + first.yz * second.yz);
}

/// @brief 从速度梯度返回旋度。
inline Vector3 curl(const Tensor3& gradient) {
    return {
        gradient[2][1] - gradient[1][2],
        gradient[0][2] - gradient[2][0],
        gradient[1][0] - gradient[0][1]};
}

} // namespace SF::Math
