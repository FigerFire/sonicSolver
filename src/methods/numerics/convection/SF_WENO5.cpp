/// @file SF_WENO5.cpp
/// @brief WENO/TENO/对流重建数值格式实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_WENO5.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace SF {
namespace WENO5 {

namespace {

void requireFiniteStencil(const double* v, int n, const char* scheme) {
    for (int s = 0; s < n; ++s) {
        if (!std::isfinite(v[s])) {
            std::cerr << "[SF FATAL] " << scheme
                      << " received non-finite scalar stencil value at index "
                      << s << ": " << v[s] << std::endl;
            std::exit(1);
        }
    }
}

double requireFiniteWeightSum(double sumA, const char* scheme) {
    if (!std::isfinite(sumA) || sumA <= 0.0) {
        std::cerr << "[SF FATAL] " << scheme
                  << " produced invalid nonlinear weight sum: "
                  << sumA << std::endl;
        std::exit(1);
    }
    return sumA;
}

} // namespace

double weno5_core(const double* v) {
    requireFiniteStencil(v, 5, "WENO5");

    double p0 = ( 2.0*v[0] - 7.0*v[1] + 11.0*v[2]) / 6.0;
    double p1 = (-1.0*v[1] + 5.0*v[2] +  2.0*v[3]) / 6.0;
    double p2 = ( 2.0*v[2] + 5.0*v[3] -  1.0*v[4]) / 6.0;

    const double d0 = 0.1, d1 = 0.6, d2 = 0.3;

    double b0 = 13.0/12.0 * std::pow(v[0] - 2.0*v[1] + v[2], 2) + 0.25 * std::pow(v[0] - 4.0*v[1] + 3.0*v[2], 2);
    double b1 = 13.0/12.0 * std::pow(v[1] - 2.0*v[2] + v[3], 2) + 0.25 * std::pow(v[1] - v[3], 2);
    double b2 = 13.0/12.0 * std::pow(v[2] - 2.0*v[3] + v[4], 2) + 0.25 * std::pow(3.0*v[2] - 4.0*v[3] + v[4], 2);

    const double eps = 1e-6;
    double a0 = d0 / std::pow(eps + b0, 2);
    double a1 = d1 / std::pow(eps + b1, 2);
    double a2 = d2 / std::pow(eps + b2, 2);
    double sumA = requireFiniteWeightSum(a0 + a1 + a2, "WENO5");

    return (a0 * p0 + a1 * p1 + a2 * p2) / sumA;
}

} // namespace WENO5
} // namespace SF
