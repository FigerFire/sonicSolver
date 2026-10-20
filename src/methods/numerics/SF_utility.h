/// @file SF_utility.h
/// @brief 可复用数值算子与格式辅助接口。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

#include "SF_field.h"
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace SF {
namespace Numerics {

inline bool finite(double x) {
    return std::isfinite(x);
}

inline double pressure(double rho, double ru, double rv, double rw, double e,
                       double gamma) {
    double invRho = 1.0 / rho;
    double ke = 0.5 * (ru * ru + rv * rv + rw * rw) * invRho;
    return (gamma - 1.0) * (e - ke);
}

inline double pressure(double rho, double ru, double rv, double rw, double e) {
    return pressure(rho, ru, rv, rw, e, 1.4);
}

inline double pressure(const Field& field, int i, int j, int k) {
    if (field.hasStateModel()) {
        return field.thermodynamicState(i, j, k).pressure;
    }
    return pressure(field(i, j, k, RHO),
                    field(i, j, k, RU),
                    field(i, j, k, RV),
                    field(i, j, k, RW),
                    field(i, j, k, E));
}

inline double requirePhysicalState(const char* context,
                                   double rho, double ru, double rv, double rw,
                                   double e, double gamma) {
    bool badConservative = !finite(rho) || !finite(ru) || !finite(rv) ||
                           !finite(rw) || !finite(e) || rho <= 0.0;
    double p = badConservative ? std::numeric_limits<double>::quiet_NaN()
                               : pressure(rho, ru, rv, rw, e, gamma);
    if (badConservative || !finite(p) || p <= 0.0) {
        std::cerr << "[SF FATAL] " << context
                  << ": invalid conservative state"
                  << " rho=" << rho
                  << " ru=" << ru
                  << " rv=" << rv
                  << " rw=" << rw
                  << " E=" << e
                  << " p=" << p
                  << std::endl;
        std::exit(1);
    }
    return p;
}

inline double requirePhysicalState(const char* context,
                                   double rho, double ru, double rv, double rw,
                                   double e) {
    return requirePhysicalState(context, rho, ru, rv, rw, e, 1.4);
}

inline double requirePhysicalState(const char* context, const Field& field, int i, int j, int k) {
    if (field.hasStateModel()) {
        try {
            const auto state = field.thermodynamicState(i, j, k);
            if (!finite(state.density) || state.density <= 0.0
                || !finite(state.pressure) || state.pressure <= 0.0
                || !finite(state.soundSpeed) || state.soundSpeed <= 0.0) {
                throw std::runtime_error("EOS closure returned an invalid state.");
            }
            return state.pressure;
        } catch (const std::exception& error) {
            std::cerr << "[SF FATAL] " << context << ": EOS state at ("
                      << i << "," << j << "," << k << ") is invalid: "
                      << error.what() << std::endl;
            std::exit(1);
        }
    }
    double rho = field(i, j, k, RHO);
    double ru = field(i, j, k, RU);
    double rv = field(i, j, k, RV);
    double rw = field(i, j, k, RW);
    double e = field(i, j, k, E);
    bool badConservative = !finite(rho) || !finite(ru) || !finite(rv) ||
                           !finite(rw) || !finite(e) || rho <= 0.0;
    double p = badConservative ? std::numeric_limits<double>::quiet_NaN()
                               : pressure(rho, ru, rv, rw, e);
    if (badConservative || !finite(p) || p <= 0.0) {
        std::cerr << "[SF FATAL] " << context
                  << ": invalid state at (" << i << "," << j << "," << k << ")"
                  << " rho=" << rho
                  << " ru=" << ru
                  << " rv=" << rv
                  << " rw=" << rw
                  << " E=" << e
                  << " p=" << p
                  << std::endl;
        std::exit(1);
    }
    return p;
}

} // namespace Numerics
} // namespace SF
