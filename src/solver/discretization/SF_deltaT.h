/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_deltaT.h
/// @brief CFL 条件时间步长估算。
///
/// 调用方式:
/// @code
///   double dt = deltaT(field, cfl);
/// @endcode

#include "SF_field.h"
#include "methods/numerics/structured/SF_structured.h"
#include "SF_utility.h"

#include <algorithm>
#include <cmath>

namespace SF {

/// 基于 CFL 条件计算局部最大容许时间步长。
/// 遍历所有内部实胞，取 min(CFL / (|U_ξ|+a_ξ + |U_η|+a_η + |U_ζ|+a_ζ))。
///
/// @param field  物理场 (含度规)
/// @param cfl    CFL 数
/// @return       全局最小 dt
inline double deltaT(const Field& field, double cfl, double gamma = 1.4) {
    double min_dt = 1e10;

    Math::forFluidInterior(field, [&](int i, int j, int k) {
        double rho = field(i, j, k, RHO);
        double ru  = field(i, j, k, RU);
        double rv  = field(i, j, k, RV);
        double rw  = field(i, j, k, RW);

        double p = field.hasEquationSet()
            ? Numerics::requirePhysicalState("deltaT", field, i, j, k)
            : Numerics::requirePhysicalState(
                "deltaT", rho, ru, rv, rw,
                field(i, j, k, E), gamma);
        double u = ru / rho;
        double v = rv / rho;
        double w = rw / rho;
        double a = std::sqrt(gamma * p / rho);

        double U_xi = std::abs(u * field.XiX(i, j, k)
                             + v * field.XiY(i, j, k)
                             + w * field.XiZ(i, j, k));
        double V_et = std::abs(u * field.EtX(i, j, k)
                             + v * field.EtY(i, j, k)
                             + w * field.EtZ(i, j, k));
        double W_ze = std::abs(u * field.ZeX(i, j, k)
                             + v * field.ZeY(i, j, k)
                             + w * field.ZeZ(i, j, k));

        double a_xi = a * Math::vecMag(field.XiX(i, j, k),
                                       field.XiY(i, j, k),
                                       field.XiZ(i, j, k));
        double a_et = a * Math::vecMag(field.EtX(i, j, k),
                                       field.EtY(i, j, k),
                                       field.EtZ(i, j, k));
        double a_ze = a * Math::vecMag(field.ZeX(i, j, k),
                                       field.ZeY(i, j, k),
                                       field.ZeZ(i, j, k));

        double sr = 0.0;
        if (Math::isDirectionActive(Math::XI))   sr += U_xi + a_xi;
        if (Math::isDirectionActive(Math::ETA))  sr += V_et + a_et;
        if (Math::isDirectionActive(Math::ZETA)) sr += W_ze + a_ze;
        double dt_local = cfl / (sr + 1e-10);
        if (dt_local < min_dt) min_dt = dt_local;
    });

    return min_dt;
}

} // namespace SF
