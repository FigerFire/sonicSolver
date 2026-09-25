/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/
/*--------------Sonic Fluid-------------------*/

#pragma once

/// @file SF_vector.h
/// @brief 通用三维向量运算接口。
///
/// 提供几何向量(Vector3)的基本代数运算:加减、数乘、点积、叉积、模长、归一化。
/// 本模块不依赖任何Field或网格结构,仅做纯数学向量计算。
/// 所有函数均为内联,可直接被math/下其他模块和求解器调用。

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

namespace SF {
namespace Math {

/// @brief 三维物理向量。
///
/// 纯数学 authority：与 Field / Mesh / MPI / Solver 完全解耦。
/// `core/state/SF_valueTypes.h` 通过 `using` 别名保持向后兼容。
struct Vector3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    Vector3() = default;
    Vector3(double xValue, double yValue, double zValue)
        : x(xValue), y(yValue), z(zValue) {}

    double& operator[](std::size_t component) {
        return component == 0 ? x : (component == 1 ? y : z);
    }
    const double& operator[](std::size_t component) const {
        return component == 0 ? x : (component == 1 ? y : z);
    }

    Vector3 operator+(const Vector3& value) const {
        return {x + value.x, y + value.y, z + value.z};
    }
    Vector3 operator-(const Vector3& value) const {
        return {x - value.x, y - value.y, z - value.z};
    }
    Vector3 operator*(double scale) const {
        return {x * scale, y * scale, z * scale};
    }
    Vector3 operator/(double scale) const {
        return {x / scale, y / scale, z / scale};
    }

    friend Vector3 operator*(double scale, const Vector3& value) {
        return value * scale;
    }
};

/// @brief 三维向量叉积。
inline Vector3 cross(const Vector3& first, const Vector3& second) {
    return {
        first.y * second.z - first.z * second.y,
        first.z * second.x - first.x * second.z,
        first.x * second.y - first.y * second.x
    };
}

/// @brief 三维向量点积。
inline double dot(const Vector3& first, const Vector3& second) {
    return first.x * second.x + first.y * second.y + first.z * second.z;
}

/// @brief 三维向量模长。
inline double norm(const Vector3& value) {
    return std::sqrt(dot(value, value));
}

/// @brief 三维向量归一化，退化时返回零向量。
inline Vector3 normalize(const Vector3& value) {
    const double magnitude = norm(value);
    return magnitude > 1.0e-12 ? value / magnitude : Vector3();
}

// ============================================================
//  向量运算自由函数
// ============================================================

/// @brief 向量点积。
/// @param a 向量a。
/// @param b 向量b。
/// @return a·b的标量结果。
/// @brief 向量归一化，退化时返回给定备用方向。
/// @param v 向量。
/// @param fallback 退化时使用的单位方向。
/// @return 单位向量；若v退化则返回fallback。
inline Vector3 normalizeOr(const Vector3& v,
                           const Vector3& fallback = Vector3(1.0, 0.0, 0.0)) {
    double n = norm(v);
    return (n > 1e-12) ? (v * (1.0 / n)) : fallback;
}

/// @brief 从C数组构造三维向量。
/// @param v 三维数组。
/// @return 三维向量。
inline Vector3 vectorFromArray(const double v[3]) {
    return Vector3(v[0], v[1], v[2]);
}

/// @brief 从std::array构造三维向量。
/// @param v 三维数组。
/// @return 三维向量。
inline Vector3 vectorFromArray(const std::array<double, 3>& v) {
    return Vector3(v[0], v[1], v[2]);
}

/// @brief 把三维向量转为std::array。
/// @param v 三维向量。
/// @return 三维数组。
inline std::array<double, 3> vectorToArray(const Vector3& v) {
    return {v.x, v.y, v.z};
}

/// @brief 原地归一化C数组形式的三维向量。
/// @param v 三维数组。
/// @param fallback 退化时使用的单位方向。
inline void normalizeArray(double v[3],
                           const Vector3& fallback = Vector3(1.0, 0.0, 0.0)) {
    const Vector3 unit = normalizeOr(vectorFromArray(v), fallback);
    v[0] = unit.x;
    v[1] = unit.y;
    v[2] = unit.z;
}

/// @brief 计算两点之间的有符号距离向量。
/// @param from 起点。
/// @param to 终点。
/// @return to - from的向量。
inline Vector3 displacement(const Vector3& from, const Vector3& to) {
    return to - from;
}

/// @brief 计算两点之间的欧几里得距离。
/// @param a 点a。
/// @param b 点b。
/// @return |b - a|的距离标量。
inline double distance(const Vector3& a, const Vector3& b) {
    return norm(b - a);
}

/// @brief 向量各分量取绝对值。
/// @param v 向量。
/// @return 分量绝对值构成的新向量。
inline Vector3 absVector(const Vector3& v) {
    return Vector3(std::abs(v.x), std::abs(v.y), std::abs(v.z));
}

/// @brief 向量按分量取最大值。
/// @param a 向量a。
/// @param b 向量b。
/// @return 逐分量取max的新向量。
inline Vector3 maxComponents(const Vector3& a, const Vector3& b) {
    return Vector3(
        std::max(a.x, b.x),
        std::max(a.y, b.y),
        std::max(a.z, b.z)
    );
}

/// @brief 向量按分量取最小值。
/// @param a 向量a。
/// @param b 向量b。
/// @return 逐分量取min的新向量。
inline Vector3 minComponents(const Vector3& a, const Vector3& b) {
    return Vector3(
        std::min(a.x, b.x),
        std::min(a.y, b.y),
        std::min(a.z, b.z)
    );
}

/// @brief 计算3×3雅可比矩阵的行列式(单元体积)。
/// @param xr ∂x/∂ξ。
/// @param yr ∂y/∂ξ。
/// @param zr ∂z/∂ξ。
/// @param xs ∂x/∂η。
/// @param ys ∂y/∂η。
/// @param zs ∂z/∂η。
/// @param xt ∂x/∂ζ。
/// @param yt ∂y/∂ζ。
/// @param zt ∂z/∂ζ。
/// @return 雅可比行列式的值J。
inline double jacobianDeterminant(
    double xr, double yr, double zr,
    double xs, double ys, double zs,
    double xt, double yt, double zt)
{
    return xr * (ys * zt - zs * yt)
         - yr * (xs * zt - zs * xt)
         + zr * (xs * yt - ys * xt);
}

/// @brief 计算度规系数ξ_x, ξ_y, ξ_z及1/J。
/// @param xr, yr, zr ∂x/∂ξ等偏导数。
/// @param xs, ys, zs ∂x/∂η等偏导数。
/// @param xt, yt, zt ∂x/∂ζ等偏导数。
/// @param invJ 1/J,供外部预计算传入。
/// @param xi_x, xi_y, xi_z 输出度规ξ分量。
/// @param et_x, et_y, et_z 输出度规η分量。
/// @param ze_x, ze_y, ze_z 输出度规ζ分量。
inline void computeMetricsFromJacobian(
    double xr, double yr, double zr,
    double xs, double ys, double zs,
    double xt, double yt, double zt,
    double invJ,
    double& xi_x, double& xi_y, double& xi_z,
    double& et_x, double& et_y, double& et_z,
    double& ze_x, double& ze_y, double& ze_z)
{
    xi_x = (ys * zt - zs * yt) * invJ;
    xi_y = (zs * xt - xs * zt) * invJ;
    xi_z = (xs * yt - ys * xt) * invJ;

    et_x = (zt * yr - yt * zr) * invJ;
    et_y = (xr * zt - zr * xt) * invJ;
    et_z = (xt * yr - xr * yt) * invJ;

    ze_x = (yr * zs - zr * ys) * invJ;
    ze_y = (zr * xs - xr * zs) * invJ;
    ze_z = (xr * ys - yr * xs) * invJ;
}

} // namespace Math
} // namespace SF
