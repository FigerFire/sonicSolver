/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_eulerState.h
/// @brief Euler守恒/原始变量转换与法向characteristic变换。
///
/// 本文件只处理局部原始变量、守恒变量和冻结法向Euler characteristic
/// variable之间的代数变换。不依赖Field、IO、MPI或求解器状态。

#include "SF_localFrame.h"

#include <array>
#include <cmath>

namespace SF {
namespace Math {
namespace Euler {

/// @brief 局部原始变量 `[rho, u_n, u_t1, u_t2, p]`。
using Primitive = std::array<double, 5>;

/// @brief 局部原始变量分量编号。
enum LocalVar {
    LRHO = 0,
    LUN  = 1,
    LUT1 = 2,
    LUT2 = 3,
    LP   = 4
};

/// @brief 从守恒变量计算全局原始变量。
/// @param q 守恒变量 `[rho,rho*u,rho*v,rho*w,E]`。
/// @param gamma 比热比。
/// @return 全局原始变量 `[rho,u,v,w,p]`。
inline Primitive primitiveFromConservative(const double q[5],
                                           double gamma = 1.4) {
    const double invRho = 1.0 / q[0];
    const double ke = 0.5 * (q[1] * q[1] + q[2] * q[2] + q[3] * q[3]) * invRho;
    return Primitive{
        q[0],
        q[1] * invRho,
        q[2] * invRho,
        q[3] * invRho,
        (gamma - 1.0) * (q[4] - ke)
    };
}

/// @brief 从全局原始变量计算守恒变量。
/// @param primitive 全局原始变量 `[rho,u,v,w,p]`。
/// @param gamma 比热比。
/// @param q 输出守恒变量 `[rho,rho*u,rho*v,rho*w,E]`。
inline void conservativeFromPrimitive(const Primitive& primitive,
                                      double gamma,
                                      double q[5]) {
    const double rho = primitive[0];
    const double u = primitive[1];
    const double v = primitive[2];
    const double w = primitive[3];
    const double p = primitive[4];

    q[0] = rho;
    q[1] = rho * u;
    q[2] = rho * v;
    q[3] = rho * w;
    q[4] = p / (gamma - 1.0) + 0.5 * rho * (u * u + v * v + w * w);
}

/// @brief 把全局原始变量转为局部原始变量。
/// @param global 全局原始变量 `[rho,u,v,w,p]`。
/// @param frame 壁面局部坐标系。
/// @return 局部原始变量 `[rho,u_n,u_t1,u_t2,p]`。
inline Primitive toLocalPrimitive(const Primitive& global,
                                  const LocalFrame& frame) {
    const Vector3 velocity(global[1], global[2], global[3]);
    return Primitive{
        global[0],
        dot(velocity, frame.n),
        dot(velocity, frame.t1),
        dot(velocity, frame.t2),
        global[4]
    };
}

/// @brief 把局部原始变量转为全局原始变量。
/// @param local 局部原始变量 `[rho,u_n,u_t1,u_t2,p]`。
/// @param frame 壁面局部坐标系。
/// @return 全局原始变量 `[rho,u,v,w,p]`。
inline Primitive toGlobalPrimitive(const Primitive& local,
                                   const LocalFrame& frame) {
    const Vector3 velocity =
        local[LUN] * frame.n + local[LUT1] * frame.t1 + local[LUT2] * frame.t2;
    return Primitive{local[LRHO], velocity.x, velocity.y, velocity.z, local[LP]};
}

/// @brief 从守恒变量直接计算局部原始变量。
/// @param q 守恒变量 `[rho,rho*u,rho*v,rho*w,E]`。
/// @param frame 壁面局部坐标系。
/// @param gamma 比热比。
/// @return 局部原始变量 `[rho,u_n,u_t1,u_t2,p]`。
inline Primitive localPrimitiveFromConservative(const double q[5],
                                                const LocalFrame& frame,
                                                double gamma = 1.4) {
    return toLocalPrimitive(primitiveFromConservative(q, gamma), frame);
}

/// @brief 计算壁面速度在局部坐标系中的分量。
/// @param wallVelocity 全局壁面速度。
/// @param frame 壁面局部坐标系。
/// @return 局部速度分量，密度和压力分量为0。
inline Primitive localWallVelocity(const double wallVelocity[3],
                                   const LocalFrame& frame) {
    const Vector3 velocity = vectorFromArray(wallVelocity);
    return Primitive{
        0.0,
        dot(velocity, frame.n),
        dot(velocity, frame.t1),
        dot(velocity, frame.t2),
        0.0
    };
}

/// @brief 由局部参考状态计算声速。
/// @param reference 局部参考原始变量。
/// @param gamma 比热比。
/// @return 参考声速；参考状态非法时返回非正数。
inline double referenceSoundSpeed(const Primitive& reference,
                                  double gamma = 1.4) {
    if (!std::isfinite(reference[LRHO]) || reference[LRHO] <= 0.0 ||
        !std::isfinite(reference[LP]) || reference[LP] <= 0.0) {
        return -1.0;
    }
    return std::sqrt(gamma * reference[LP] / reference[LRHO]);
}

/// @brief 把局部原始变量扰动投影到法向Euler characteristic variables。
/// @param value 当前样本的局部原始变量。
/// @param reference 冻结变换所用的局部壁面参考状态。
/// @param gamma 比热比。
/// @param characteristic 输出 characteristic-variable 扰动。
/// @return 参考状态物理有效时返回true。
inline bool primitiveToCharacteristicDifference(const Primitive& value,
                                                const Primitive& reference,
                                                double gamma,
                                                Primitive& characteristic) {
    const double c = referenceSoundSpeed(reference, gamma);
    if (!(c > 0.0) || !std::isfinite(c)) return false;

    const double rho = reference[LRHO];
    const double c2 = c * c;
    const double drho = value[LRHO] - reference[LRHO];
    const double dun = value[LUN] - reference[LUN];
    const double dut1 = value[LUT1] - reference[LUT1];
    const double dut2 = value[LUT2] - reference[LUT2];
    const double dp = value[LP] - reference[LP];

    characteristic[LRHO] = dp - c2 * drho;
    characteristic[LUN] = dp + rho * c * dun;
    characteristic[LUT1] = dut1;
    characteristic[LUT2] = dut2;
    characteristic[LP] = dp - rho * c * dun;
    return true;
}

/// @brief 把法向 characteristic-variable 导数还原为局部原始变量导数。
/// @param characteristicDerivative characteristic-variable 法向导数。
/// @param reference 冻结变换所用的局部壁面参考状态。
/// @param gamma 比热比。
/// @param primitiveDerivative 输出局部原始变量法向导数。
/// @return 参考状态物理有效时返回true。
inline bool characteristicDerivativeToPrimitive(
        const Primitive& characteristicDerivative,
        const Primitive& reference,
        double gamma,
        Primitive& primitiveDerivative) {
    const double c = referenceSoundSpeed(reference, gamma);
    if (!(c > 0.0) || !std::isfinite(c)) return false;

    const double rho = reference[LRHO];
    const double c2 = c * c;
    const double dp = 0.5 * (characteristicDerivative[LUN]
                           + characteristicDerivative[LP]);
    const double dun = (characteristicDerivative[LUN]
                      - characteristicDerivative[LP])
                     / (2.0 * rho * c);

    primitiveDerivative[LP] = dp;
    primitiveDerivative[LUN] = dun;
    primitiveDerivative[LRHO] = (dp - characteristicDerivative[LRHO]) / c2;
    primitiveDerivative[LUT1] = characteristicDerivative[LUT1];
    primitiveDerivative[LUT2] = characteristicDerivative[LUT2];
    return true;
}

} // namespace Euler
} // namespace Math
} // namespace SF
