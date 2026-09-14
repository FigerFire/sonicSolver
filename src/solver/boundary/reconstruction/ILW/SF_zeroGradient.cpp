/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

/// @file SF_zeroGradient.cpp
/// @brief 3D 高阶 ILW 零梯度边界条件。

#include "SF_zeroGradient.h"

#include "SF_boundaryClosure.h"

#include <array>

namespace SF {
namespace Boundary {
namespace ILW {
namespace ZeroGradient {
namespace {

using BoundaryClosure::BoundaryNormal;

void requireNormal(const BoundaryNormal& normal,
                   int i, int j, int k,
                   const char* bcName) {
    if (normal.sign == 0 || normal.axis < 0) {
        BoundaryClosure::fatalClosure(
            bcName, i, j, k,
            "point is not on an active structured physical boundary.");
    }
}

template <typename Getter, typename Setter>
void writeZeroGradientGhosts(Field& field,
                             int i, int j, int k,
                             int accuracyOrder,
                             Getter&& getter,
                             Setter&& setter,
                             const char* bcName) {
    const BoundaryNormal normal = BoundaryClosure::boundaryNormal(field, i, j, k);
    requireNormal(normal, i, j, k, bcName);

    /// 仅当该点在活跃法向上是虚胞时才做 Taylor 外推。
    const bool isNormalGhost = [&]() -> bool {
        const int ng = field.NG();
        if (normal.axis == 0)
            return i < ng || i >= field.NX() + ng;
        if (normal.axis == 1)
            return j < ng || j >= field.NY() + ng;
        if (normal.axis == 2)
            return k < ng || k >= field.NZ() + ng;
        return false;
    }();

    if (isNormalGhost) {
        int ri = 0, rj = 0, rk = 0;
        nearestInteriorCell(field, i, j, k, ri, rj, rk);

        std::array<double, BoundaryClosure::kMaxAccuracyOrder> coeff{};
        BoundaryClosure::buildMultiDimScalarTaylorCoefficients(
            field, ri, rj, rk, normal, accuracyOrder,
            getter, bcName, coeff);
        coeff[1] = 0.0;

        const int layer = BoundaryClosure::ghostLayer(i, j, k, ri, rj, rk, normal);
        if (layer < 1 || layer > field.NG()) {
            BoundaryClosure::fatalClosure(bcName, i, j, k,
                                  "invalid ghost layer for zero-gradient closure.");
        }
        const double h = BoundaryClosure::normalSpacing(field, ri, rj, rk, normal);
        setter(i, j, k,
               BoundaryClosure::evaluateTaylor(
                   coeff, BoundaryClosure::taylorOrder(accuracyOrder),
                   BoundaryClosure::ghostDistance(layer, h)));
        return;
    }

    /// 实边界点：保持原值，并向该法向的虚胞层广播。
    auto broadcast = [&](int gi, int gj, int gk,
                         int ri, int rj, int rk) {
        int refI = ri, refJ = rj, refK = rk;
        nearestInteriorCell(field, ri, rj, rk, refI, refJ, refK);

        std::array<double, BoundaryClosure::kMaxAccuracyOrder> coeff{};
        BoundaryClosure::buildMultiDimScalarTaylorCoefficients(
            field, refI, refJ, refK, normal, accuracyOrder,
            getter, bcName, coeff);
        coeff[1] = 0.0;

        const int layer =
            BoundaryClosure::ghostLayer(gi, gj, gk, refI, refJ, refK, normal);
        if (layer < 1 || layer > field.NG()) {
            BoundaryClosure::fatalClosure(bcName, gi, gj, gk,
                                  "invalid ghost layer for zero-gradient closure.");
        }

        const double h = BoundaryClosure::normalSpacing(field, refI, refJ, refK, normal);
        setter(gi, gj, gk,
               BoundaryClosure::evaluateTaylor(
                   coeff, BoundaryClosure::taylorOrder(accuracyOrder),
                   BoundaryClosure::ghostDistance(layer, h)));
    };
    forBoundaryGhostsAlongAxis(field, i, j, k, normal.axis, broadcast);
}

} // namespace

void applyScalar(Field& field, int i, int j, int k, int vIdx,
                 int accuracyOrder) {
    const BoundaryClosure::FieldVariableGetter getter{vIdx};
    auto setter = [&](int gi, int gj, int gk, double value) {
        field(gi, gj, gk, vIdx) = value;
    };
    writeZeroGradientGhosts(field, i, j, k, accuracyOrder,
                            getter, setter, "zeroGradient");
}

void applyVector3(Field& field, int i, int j, int k, int vIdx,
                  int accuracyOrder) {
    for (int comp = 0; comp < 3; ++comp) {
        const BoundaryClosure::FieldVariableGetter getter{vIdx + comp};
        auto setter = [&](int gi, int gj, int gk, double value) {
            field(gi, gj, gk, vIdx + comp) = value;
        };
        writeZeroGradientGhosts(field, i, j, k, accuracyOrder,
                                getter, setter,
                                "zeroGradientVector3");
    }
}

} // namespace ZeroGradient
} // namespace ILW
} // namespace Boundary
} // namespace SF
