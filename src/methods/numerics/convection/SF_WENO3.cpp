/// @file SF_WENO3.cpp
/// @brief WENO/TENO/对流重建数值格式实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_WENO3.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace SF {
namespace WENO3 {

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

double weno3_core(const double* v) {
    requireFiniteStencil(v, 4, "WENO3");

    double p0 = (-v[0] + 5.0*v[1] + 2.0*v[2]) / 6.0;
    double p1 = ( 2.0*v[1] + 5.0*v[2] -   v[3]) / 6.0;

    const double d0 = 1.0/3.0, d1 = 2.0/3.0;
    double b0 = (v[2] - v[0]) * (v[2] - v[0]);
    double b1 = (v[3] - v[1]) * (v[3] - v[1]);

    const double eps = 1e-6;
    double a0 = d0 / ((eps + b0) * (eps + b0));
    double a1 = d1 / ((eps + b1) * (eps + b1));
    double sumA = requireFiniteWeightSum(a0 + a1, "WENO3");
    return (a0 * p0 + a1 * p1) / sumA;
}

} // namespace WENO3
} // namespace SF
