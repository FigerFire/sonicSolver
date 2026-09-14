/// @file SF_fixedValue.cpp
/// @brief 线性 ghost 层的 fixedValue 标量/向量实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_fixedValue.h"
#include "core/mesh/SF_dimension.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace SF {
namespace Boundary {
namespace Algebraic {
namespace FixedValue {
namespace {

bool isGhostCell(const Field& field, int i, int j, int k) {
    const int ng = field.NG();
    return (i < ng || i >= field.NX() + ng ||
            j < ng || j >= field.NY() + ng ||
            k < ng || k >= field.NZ() + ng);
}

void nearestInteriorCell(const Field& field,
                         int i, int j, int k,
                         int& ri, int& rj, int& rk) {
    const int ng = field.NG();
    ri = std::max(ng, std::min(i, field.NX() + ng - 1));
    rj = std::max(ng, std::min(j, field.NY() + ng - 1));
    rk = std::max(ng, std::min(k, field.NZ() + ng - 1));
}

template <typename Fn>
void forBoundaryGhostsAlongAxis(Field& field,
                                int i, int j, int k,
                                int axis,
                                Fn&& fn) {
    const int ng = field.NG();
    if (!Math::isDirectionActiveIndex(axis)) return;

    if (axis == 0) {
        if (i == ng) {
            for (int b = 0; b < ng; ++b) fn(b, j, k, i, j, k);
        }
        if (i == field.NX() + ng - 1) {
            for (int b = 1; b <= ng; ++b) fn(i + b, j, k, i, j, k);
        }
    } else if (axis == 1) {
        if (j == ng) {
            for (int b = 0; b < ng; ++b) fn(i, b, k, i, j, k);
        }
        if (j == field.NY() + ng - 1) {
            for (int b = 1; b <= ng; ++b) fn(i, j + b, k, i, j, k);
        }
    } else if (axis == 2) {
        if (k == ng) {
            for (int b = 0; b < ng; ++b) fn(i, j, b, i, j, k);
        }
        if (k == field.NZ() + ng - 1) {
            for (int b = 1; b <= ng; ++b) fn(i, j, k + b, i, j, k);
        }
    }
}

double requireDensity(const Field& field,
                      int i, int j, int k,
                      const char* context) {
    const double rho = field.hasEquationSet()
        ? field.thermodynamicState(i, j, k).density
        : field(i, j, k, RHO);
    if (!std::isfinite(rho) || rho <= 0.0) {
        std::ostringstream oss;
        oss << context << ": invalid rho at (" << i << "," << j << ","
            << k << "), rho=" << rho;
        throw std::runtime_error(oss.str());
    }
    return rho;
}

SF::Vector3 primitiveVelocity(const Field& field,
                              int i, int j, int k,
                              const char* context) {
    const double rho = requireDensity(field, i, j, k, context);
    const int momentum = field.hasEquationSet()
        ? field.equationSet()->momentumIndex(0) : RU;
    return SF::Vector3(field(i, j, k, momentum) / rho,
                       field(i, j, k, momentum+1) / rho,
                       field(i, j, k, momentum+2) / rho);
}

void setMomentumFromVelocity(Field& field,
                             int i, int j, int k,
    const SF::Vector3& velocity,
                             const char* context) {
    const double rho = requireDensity(field, i, j, k, context);
    const int momentum = field.hasEquationSet()
        ? field.equationSet()->momentumIndex(0) : RU;
    field(i, j, k, momentum) = rho * velocity.x;
    field(i, j, k, momentum+1) = rho * velocity.y;
    field(i, j, k, momentum+2) = rho * velocity.z;
}

} // namespace

void applyScalar(Field& field, int i, int j, int k, int axis,
                 double bcValue, int vIdx) {
    auto applyFixed = [&](int gi, int gj, int gk,
                          int ri, int rj, int rk) {
        field(gi, gj, gk, vIdx) =
            2.0 * bcValue - field(ri, rj, rk, vIdx);
    };

    if (isGhostCell(field, i, j, k)) {
        int ri = 0, rj = 0, rk = 0;
        nearestInteriorCell(field, i, j, k, ri, rj, rk);
        applyFixed(i, j, k, ri, rj, rk);
        return;
    }

    field(i, j, k, vIdx) = bcValue;
    forBoundaryGhostsAlongAxis(field, i, j, k, axis, applyFixed);
}

void applyVector3(Field& field, int i, int j, int k, int axis,
                  const SF::Vector3& bcValue, int vIdx) {
    const int momentum = field.hasEquationSet()
        ? field.equationSet()->momentumIndex(0) : RU;
    if (vIdx != momentum) {
        throw std::runtime_error(
            "fixedValue Vector3 boundary is only implemented for velocity/RU.");
    }
    auto applyFixed = [&](int gi, int gj, int gk,
                          int ri, int rj, int rk) {
        const SF::Vector3 realValue =
            primitiveVelocity(field, ri, rj, rk,
                              "fixedValue velocity boundary reference");
        setMomentumFromVelocity(
            field, gi, gj, gk, 2.0 * bcValue - realValue,
            "fixedValue velocity boundary ghost");
    };

    if (isGhostCell(field, i, j, k)) {
        int ri = 0, rj = 0, rk = 0;
        nearestInteriorCell(field, i, j, k, ri, rj, rk);
        applyFixed(i, j, k, ri, rj, rk);
        return;
    }

    setMomentumFromVelocity(field, i, j, k, bcValue,
                            "fixedValue velocity boundary");
    forBoundaryGhostsAlongAxis(field, i, j, k, axis, applyFixed);
}

} // namespace FixedValue
} // namespace Algebraic
} // namespace Boundary
} // namespace SF
