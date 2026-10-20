/// @file SF_WENO7.cpp
/// @brief WENO/TENO/对流重建数值格式实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_WENO7.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace SF {
namespace WENO7 {

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

// ── 辅助: 4-point stencil 平滑度指示器 ──
static double beta4(const double v[4], double a, double b) {
    double d1 = v[1] - v[0];
    double d2 = (v[2] - 2.0*v[1] + v[0]) * 0.5;
    double d3 = (v[3] - 3.0*v[2] + 3.0*v[1] - v[0]) / 6.0;

    double A = d1 - d2 + 2.0*d3;
    double B = 2.0*d2 - 6.0*d3;
    double C = 3.0*d3;

    auto I1 = [&](auto f1) { return f1(b) - f1(a); };

    double beta1 = A*A * I1([](double x){ return x; })
                 + A*B * I1([](double x){ return x*x; })
                 + (B*B + 2.0*A*C) * I1([](double x){ return x*x*x/3.0; })
                 + B*C * I1([](double x){ double x2=x*x; return x2*x2/2.0; })
                 + C*C * I1([](double x){ return x*x*x*x*x/5.0; });

    double D = 2.0*d2 - 6.0*d3;
    double E = 6.0*d3;
    double beta2 = D*D * I1([](double x){ return x; })
                 + D*E * I1([](double x){ return x*x; })
                 + E*E * I1([](double x){ return x*x*x/3.0; });

    double beta3 = 36.0 * d3 * d3 * I1([](double x){ return x; });
    return beta1 + beta2 + beta3;
}

double weno7_core(const double* v) {
    requireFiniteStencil(v, 7, "WENO7");

    double p0 = (-3.0*v[0] + 13.0*v[1] - 23.0*v[2] + 25.0*v[3]) / 12.0;
    double p1 = (       v[1] -  5.0*v[2] + 13.0*v[3] +  3.0*v[4]) / 12.0;
    double p2 = (     - v[2] +  7.0*v[3] +  7.0*v[4] -      v[5]) / 12.0;
    double p3 = ( 3.0*v[3] + 13.0*v[4] -  5.0*v[5] +      v[6]) / 12.0;

    const double d0 = 1.0/35.0, d1 = 12.0/35.0, d2 = 18.0/35.0, d3 = 4.0/35.0;

    double sub0[4] = {v[0], v[1], v[2], v[3]};
    double sub1[4] = {v[1], v[2], v[3], v[4]};
    double sub2[4] = {v[2], v[3], v[4], v[5]};
    double sub3[4] = {v[3], v[4], v[5], v[6]};

    double b0 = beta4(sub0, 2.5, 3.5);
    double b1 = beta4(sub1, 1.5, 2.5);
    double b2 = beta4(sub2, 0.5, 1.5);
    double b3 = beta4(sub3, -0.5, 0.5);

    const double eps = 1e-6;
    double a0 = d0 / ((eps + b0) * (eps + b0));
    double a1 = d1 / ((eps + b1) * (eps + b1));
    double a2 = d2 / ((eps + b2) * (eps + b2));
    double a3 = d3 / ((eps + b3) * (eps + b3));
    double sumA = requireFiniteWeightSum(a0 + a1 + a2 + a3, "WENO7");

    return (a0*p0 + a1*p1 + a2*p2 + a3*p3) / sumA;
}

} // namespace WENO7
} // namespace SF
