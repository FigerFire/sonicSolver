/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_localFrame.h
/// @brief 局部正交坐标系与坐标投影工具。
///
/// 本文件只处理三维局部坐标系构造、点到局部坐标的投影和壁面有符号距离。
/// 不依赖Field、IO、MPI或求解器状态。

#include "SF_vector.h"

#include <array>
#include <cmath>

namespace SF {
namespace Math {

/// @brief 壁面局部正交坐标系。
struct LocalFrame {
    /// @brief 局部法向。
    Vector3 n = Vector3(1.0, 0.0, 0.0);
    /// @brief 第一切向。
    Vector3 t1 = Vector3(0.0, 1.0, 0.0);
    /// @brief 第二切向。
    Vector3 t2 = Vector3(0.0, 0.0, 1.0);
};

/// @brief 构造壁面局部正交坐标系。
/// @param wallNormal 局部法向。
/// @return 由法向和两个切向组成的局部坐标系。
inline LocalFrame makeLocalFrame(const double wallNormal[3]) {
    LocalFrame frame;
    frame.n = normalizeOr(vectorFromArray(wallNormal));

    Vector3 ref(1.0, 0.0, 0.0);
    if (std::abs(dot(ref, frame.n)) > 0.9) {
        ref = Vector3(0.0, 1.0, 0.0);
    }

    const double projection = dot(ref, frame.n);
    frame.t1 = normalizeOr(ref - projection * frame.n, Vector3(0.0, 1.0, 0.0));
    frame.t2 = normalizeOr(cross(frame.n, frame.t1), Vector3(0.0, 0.0, 1.0));
    return frame;
}

/// @brief 构造二维局部坐标系。
///
/// inactive方向被固定为第二切向，壁面法向先投影到活动平面，
/// 再用活动平面内的第一切向构成局部二维 `(s,a)` 拟合坐标。
///
/// @param wallNormal 壁面局部法向。
/// @param inactiveNormal 关闭方向在物理空间中的单位法向。
/// @return 适用于二维局部拟合的坐标系。
inline LocalFrame makeLocalFrame2D(const double wallNormal[3],
                                   const std::array<double, 3>& inactiveNormal) {
    const Vector3 t2 = normalizeOr(vectorFromArray(inactiveNormal),
                                   Vector3(0.0, 0.0, 1.0));
    const Vector3 normal = vectorFromArray(wallNormal);
    const double normalProjection = dot(normal, t2);
    const Vector3 projected = normal - normalProjection * t2;
    const double mag2 = dot(projected, projected);
    if (mag2 <= 1e-28) {
        return makeLocalFrame(wallNormal);
    }

    LocalFrame frame;
    frame.n = normalizeOr(projected);
    frame.t2 = t2;
    frame.t1 = normalizeOr(cross(frame.t2, frame.n), Vector3(0.0, 1.0, 0.0));
    frame.t2 = normalizeOr(cross(frame.n, frame.t1), Vector3(0.0, 0.0, 1.0));
    return frame;
}

/// @brief 计算点到壁面的有符号法向距离。
/// @param p 空间点。
/// @param wallPoint 壁面参考点。
/// @param wallNormal 壁面局部法向。
/// @return `(p-wallPoint) dot wallNormal`。
inline double signedDistanceFromWall(const double p[3],
                                     const double wallPoint[3],
                                     const double wallNormal[3]) {
    return dot(vectorFromArray(p) - vectorFromArray(wallPoint),
               vectorFromArray(wallNormal));
}

/// @brief 计算点在壁面局部坐标系中的坐标。
/// @param p 空间点。
/// @param wallPoint 壁面参考点。
/// @param frame 壁面局部坐标系。
/// @param s 输出法向坐标。
/// @param a 输出第一切向坐标。
/// @param b 输出第二切向坐标。
inline void localCoordinates(const double p[3],
                             const double wallPoint[3],
                             const LocalFrame& frame,
                             double& s,
                             double& a,
                             double& b) {
    const Vector3 r = vectorFromArray(p) - vectorFromArray(wallPoint);
    s = dot(r, frame.n);
    a = dot(r, frame.t1);
    b = dot(r, frame.t2);
}

} // namespace Math
} // namespace SF
