#pragma once

/// @file SF_vectorCalculus.h
/// @brief 结构曲线网格上的共享梯度、速度梯度和旋度算子。

#include "core/mesh/SF_dimension.h"
#include "SF_field.h"
#include "SF_scalarField.h"
#include "methods/math/SF_tensor.h"

#include <array>

namespace SF::Numerics::VectorCalculus {

inline void offset(
        int axis, int sign, int& di, int& dj, int& dk) {
    di = dj = dk = 0;
    if (axis == 0) di = sign;
    else if (axis == 1) dj = sign;
    else dk = sign;
}

inline std::array<double, 3> metric(
        const Field& field, int axis, int i, int j, int k) {
    if (axis == 0) {
        return {field.XiX(i,j,k), field.XiY(i,j,k), field.XiZ(i,j,k)};
    }
    if (axis == 1) {
        return {field.EtX(i,j,k), field.EtY(i,j,k), field.EtZ(i,j,k)};
    }
    return {field.ZeX(i,j,k), field.ZeY(i,j,k), field.ZeZ(i,j,k)};
}

/// @brief 把计算空间三个方向导数转换为物理空间梯度。
inline Vector3 physicalGradient(
        const Field& field,
        double dXi,
        double dEta,
        double dZeta,
        int i,
        int j,
        int k) {
    return {
        dXi * field.XiX(i,j,k)
            + dEta * field.EtX(i,j,k)
            + dZeta * field.ZeX(i,j,k),
        dXi * field.XiY(i,j,k)
            + dEta * field.EtY(i,j,k)
            + dZeta * field.ZeY(i,j,k),
        dXi * field.XiZ(i,j,k)
            + dEta * field.EtZ(i,j,k)
            + dZeta * field.ZeZ(i,j,k)};
}

} // namespace SF::Numerics::VectorCalculus

namespace SF::CENTRAL2 {

/// @brief 二阶中心梯度需要的一侧 halo 深度。
inline constexpr int requiredHaloDepth = 1;

/// @brief 对任意标量采样器计算二阶中心物理梯度。
/// @tparam Sample 可调用对象，签名为 `double(int,int,int)`。
/// @note 调用前必须由 Discretization 完成边界闭合和 halo 深度保证。
template <typename Sample>
inline Vector3 gradSampled(
        const Field& field,
        Sample&& sample,
        int i,
        int j,
        int k) {
    double derivatives[3]{0.0, 0.0, 0.0};
    for (int axis = 0; axis < 3; ++axis) {
        if (!Math::isDirectionActiveIndex(axis)) continue;
        int ip = 0, jp = 0, kp = 0;
        int im = 0, jm = 0, km = 0;
        Numerics::VectorCalculus::offset(axis, 1, ip, jp, kp);
        Numerics::VectorCalculus::offset(axis, -1, im, jm, km);
        derivatives[axis] = 0.5 * (
            sample(i + ip, j + jp, k + kp)
            - sample(i + im, j + jm, k + km));
    }
    return Numerics::VectorCalculus::physicalGradient(
        field, derivatives[0], derivatives[1], derivatives[2], i, j, k);
}

/// @brief 二阶中心差分标量梯度。
inline Vector3 grad(
        const Field& field,
        const ScalarField& scalar,
        int i,
        int j,
        int k) {
    return gradSampled(
        field,
        [&](int ii, int jj, int kk) { return scalar(ii, jj, kk); },
        i, j, k);
}

/// @brief 二阶中心差分向量梯度，`result[a][b]=dU_a/dx_b`。
inline Tensor3 grad(
        const Field& field,
        const std::array<ScalarField, 3>& vector,
        int i,
        int j,
        int k) {
    Tensor3 result;
    for (int component = 0; component < 3; ++component) {
        result[(size_t)component] = grad(
            field, vector[(size_t)component], i, j, k);
    }
    return result;
}

/// @brief 二阶中心差分速度旋度 `curl(U)`。
inline Vector3 curl(
        const Field& field,
        const std::array<ScalarField, 3>& velocity,
        int i,
        int j,
        int k) {
    return Math::curl(grad(field, velocity, i, j, k));
}

} // namespace SF::CENTRAL2
