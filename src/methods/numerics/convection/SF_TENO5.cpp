/// @file SF_TENO5.cpp
/// @brief WENO/TENO/对流重建数值格式实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_TENO5.h"

#include <atomic>
#include <cmath>
#include <cstdlib>
#include <iostream>

namespace SF {
namespace TENO5 {

namespace {

std::atomic<unsigned long long> teno5AllCutoffCount{0};

void requireFiniteStencil(const double* v, int n) {
    for (int s = 0; s < n; ++s) {
        if (!std::isfinite(v[s])) {
            std::cerr << "[SF FATAL] TENO5 received non-finite scalar stencil "
                      << "value at index " << s << ": " << v[s]
                      << std::endl;
            std::exit(1);
        }
    }
}

double requireFiniteGammaSum(double gammaSum) {
    if (!std::isfinite(gammaSum) || gammaSum <= 0.0) {
        std::cerr << "[SF FATAL] TENO5 produced invalid gamma sum: "
                  << gammaSum << std::endl;
        std::exit(1);
    }
    return gammaSum;
}

} // namespace

double teno5_core(const double* v) {
    requireFiniteStencil(v, 5);

    const double p0 = ( 2.0 * v[0] - 7.0 * v[1] + 11.0 * v[2]) / 6.0;
    const double p1 = (-1.0 * v[1] + 5.0 * v[2] +  2.0 * v[3]) / 6.0;
    const double p2 = ( 2.0 * v[2] + 5.0 * v[3] -  1.0 * v[4]) / 6.0;

    const double b0 =
        13.0 / 12.0 * std::pow(v[0] - 2.0 * v[1] + v[2], 2)
        + 0.25 * std::pow(v[0] - 4.0 * v[1] + 3.0 * v[2], 2);
    const double b1 =
        13.0 / 12.0 * std::pow(v[1] - 2.0 * v[2] + v[3], 2)
        + 0.25 * std::pow(v[1] - v[3], 2);
    const double b2 =
        13.0 / 12.0 * std::pow(v[2] - 2.0 * v[3] + v[4], 2)
        + 0.25 * std::pow(3.0 * v[2] - 4.0 * v[3] + v[4], 2);

    constexpr double d0 = 0.1;
    constexpr double d1 = 0.6;
    constexpr double d2 = 0.3;
    constexpr double eps = 1.0e-40;
    constexpr double ct = 1.0e-5;
    constexpr int q = 6;

    const double tau5 = std::abs(b0 - b2);
    const double gamma0 = std::pow(1.0 + tau5 / (b0 + eps), q);
    const double gamma1 = std::pow(1.0 + tau5 / (b1 + eps), q);
    const double gamma2 = std::pow(1.0 + tau5 / (b2 + eps), q);
    const double gammaSum = requireFiniteGammaSum(gamma0 + gamma1 + gamma2);

    const double c0 = (gamma0 / gammaSum < ct) ? 0.0 : 1.0;
    const double c1 = (gamma1 / gammaSum < ct) ? 0.0 : 1.0;
    const double c2 = (gamma2 / gammaSum < ct) ? 0.0 : 1.0;
    const double activeWeight = d0 * c0 + d1 * c1 + d2 * c2;

    if (activeWeight <= 0.0) {
        const unsigned long long count =
            ++teno5AllCutoffCount;
        std::cerr << "[SF FATAL] TENO5 rejected all substencils; "
                  << "all-cutoff count=" << count
                  << ". This would silently revert to a linear stencil. "
                  << "Use WENO5 for this case or adjust TENO cutoff "
                  << "parameters after checking the local discontinuity."
                  << std::endl;
        std::exit(1);
    }
    return (d0 * c0 * p0 + d1 * c1 * p1 + d2 * c2 * p2)
         / activeWeight;
}

} // namespace TENO5
} // namespace SF
