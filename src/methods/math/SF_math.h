/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/
/*--------------Sonic Fluid-------------------*/

#pragma once

/// @file SF_math.h
/// @brief 通用数学模块统一接口。
///
/// 本文件聚合math/下所有纯数学操作, 供求解器、IBM、边界条件等模块调用。
/// 所有函数与Field/网格拓扑解耦, 仅依赖数学数据结构。
///
/// 包含:
/// - SF_vector.h   : 三维向量(Vector3)代数运算
/// - SF_tensor.h   : 三维张量、对称张量和缩并运算
/// - SF_taylor.h   : Taylor展开与外推
/// - SF_derivative.h: 非均匀有限差分权重与导数计算
/// - SF_polynomial.h: 局部多项式基与最小二乘投影矩阵
/// - SF_linearSystem.h: 小型稠密线性方程组求解
/// - SF_localFrame.h: 局部坐标系与坐标投影
/// - SF_surfaceCurvature.h: 局部曲率张量拟合
/// - SF_weightedProjection.h: 向量样本权重代入与投影指标
///
/// 使用方式:
/// @code
///   #include "SF_math.h"
///   using namespace SF::Math;
///
///   Vector3 v = cross(a, b);        // 叉积
///   double d  = dot(a, b);           // 点积
///   auto q    = Taylor::evaluate1D<5>(derivs, dist, order); // Taylor外推
/// @endcode

#include "SF_vector.h"
#include "SF_tensor.h"
#include "SF_taylor.h"
#include "SF_derivative.h"
#include "methods/math/discrete/SF_polynomial.h"
#include "SF_linearSystem.h"
#include "discrete/SF_denseLU.h"
#include "SF_localFrame.h"
#include "SF_surfaceCurvature.h"
#include "SF_weightedProjection.h"

namespace SF {
namespace Math {

/// @brief 一阶前向差分。
/// @param qp 正侧值 q_{i+1}。
/// @param qc 当前值 q_i。
/// @param h 网格间距。
/// @return (q_{i+1} - q_i) / h。
inline double forwardDiff(double qp, double qc, double h) {
    return (qp - qc) / h;
}

/// @brief 一阶后向差分。
/// @param qc 当前值 q_i。
/// @param qm 负侧值 q_{i-1}。
/// @param h 网格间距。
/// @return (q_i - q_{i-1}) / h。
inline double backwardDiff(double qc, double qm, double h) {
    return (qc - qm) / h;
}

/// @brief 二阶中心差分。
/// @param qp 正侧值 q_{i+1}。
/// @param qm 负侧值 q_{i-1}。
/// @param h 网格间距。
/// @return (q_{i+1} - q_{i-1}) / (2h)。
inline double centralDiff(double qp, double qm, double h) {
    return (qp - qm) / (2.0 * h);
}

/// @brief 四点一阶非均匀差分(用于边界附近)。
/// @param q0 最近点值。
/// @param q1 第二近点值。
/// @param q2 第三近点值。
/// @param h0, h1, h2 依次间距。
/// @return 一阶导数的非均匀估计。
inline double nonuniformFirstDiff(
    double q0, double q1, double q2,
    double h0, double h1, double h2)
{
    double h01 = h0 + h1;
    double h02 = h0 + h2;
    return ((h01 + h02) / (h01 * h02) * q1
          - h02 / (h01 * (h02 - h01)) * q2
          + (h01 + h02 - 2.0 * h0) / (h01 * h02) * q0);
}

/// @brief 限制器: minmod函数。
/// @param a 第一个参数。
/// @param b 第二个参数。
/// @return sign(a)*min(|a|,|b|) 如果a和b同号; 否则返回0。
inline double minmod(double a, double b) {
    if (a * b <= 0.0) return 0.0;
    return (std::abs(a) < std::abs(b)) ? a : b;
}

/// @brief 速度平方 (RU²+RV²+RW²) / RHO²。
/// @param rho 密度。
/// @param ru x方向动量。
/// @param rv y方向动量。
/// @param rw z方向动量。
/// @return 速度模长的平方 u²+v²+w²。
inline double velocitySquared(double rho, double ru, double rv, double rw) {
    double invRho = 1.0 / rho;
    return (ru * ru + rv * rv + rw * rw) * invRho * invRho;
}

/// @brief 从守恒量计算压力 p = (γ-1)(E - KE)。
/// @param rho 密度。
/// @param ru x方向动量。
/// @param rv y方向动量。
/// @param rw z方向动量。
/// @param e 总能量。
/// @param gamma 比热比, 默认1.4。
/// @return 压力 p。
inline double pressureFromConserved(
    double rho, double ru, double rv, double rw, double e,
    double gamma = 1.4)
{
    double invRho = 1.0 / rho;
    double ke = 0.5 * (ru * ru + rv * rv + rw * rw) * invRho;
    return (gamma - 1.0) * (e - ke);
}

/// @brief 从守恒量计算声速 c = sqrt(γ p / ρ)。
/// @param rho 密度。
/// @param ru x方向动量。
/// @param rv y方向动量。
/// @param rw z方向动量。
/// @param e 总能量。
/// @param gamma 比热比, 默认1.4。
/// @return 声速 c。
inline double soundSpeed(
    double rho, double ru, double rv, double rw, double e,
    double gamma = 1.4)
{
    double p = pressureFromConserved(rho, ru, rv, rw, e, gamma);
    return std::sqrt(gamma * p / rho);
}

/// @brief 判断值是否有限(非NaN非Inf)。
/// @param x 待检值。
/// @return true表示有限。
inline bool isFinite(double x) {
    return std::isfinite(x);
}

} // namespace Math
} // namespace SF
