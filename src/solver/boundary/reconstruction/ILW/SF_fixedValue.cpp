/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

/// @file SF_fixedValue.cpp
/// @brief 3D 高阶 ILW 固定值边界条件。
///
/// 对结构网格边界使用多维张量积 WENO（`buildMultiDimScalarTaylorCoefficients`）
/// 计算法向 Taylor 系数，样本不足或拟合失败时 fail-fast。
/// 固定值 BC 覆盖系数常数项为预设值，然后 Taylor 展开到 ghost 层。

#include "SF_fixedValue.h"

#include "SF_boundaryClosure.h"

#include <array>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace SF {
namespace Boundary {
namespace ILW {
namespace FixedValue {
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

double requireDensity(const Field& field,
                      int i, int j, int k,
                      const char* context) {
    const double rho = field(i, j, k, RHO);
    if (!std::isfinite(rho) || rho <= 0.0) {
        std::ostringstream oss;
        oss << context << ": invalid rho at (" << i << "," << j << ","
            << k << "), rho=" << rho;
        throw std::runtime_error(oss.str());
    }
    return rho;
}

struct PrimitiveVelocityComponentGetter {
    int comp = 0;

    double operator()(const Field& field, int i, int j, int k) const {
        const double rho =
            requireDensity(field, i, j, k,
                           "ILW fixedValue velocity sample");
        return field(i, j, k, RU + comp) / rho;
    }
};

template <typename Setter, typename Getter>
void writeFixedValueGhosts(Field& field,
                           int i, int j, int k,
                           double bcValue,
                           int accuracyOrder,
                           const char* bcName,
                           Setter&& setter,
                           const Getter& getter) {
    const BoundaryNormal normal = BoundaryClosure::boundaryNormal(field, i, j, k);
    requireNormal(normal, i, j, k, bcName);

    /// 仅当该点在活跃法向上是虚胞时才做 Taylor 外推；
    /// 2D 问题中 empty 方向的虚胞会由 EMPTY BC 另行处理。
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
        coeff[0] = bcValue;

        const int layer = BoundaryClosure::ghostLayer(i, j, k, ri, rj, rk, normal);
        if (layer < 1 || layer > field.NG()) {
            BoundaryClosure::fatalClosure(bcName, i, j, k,
                                  "invalid ghost layer for fixed-value closure.");
        }
        const double h = BoundaryClosure::normalSpacing(field, ri, rj, rk, normal);
        setter(i, j, k,
               BoundaryClosure::evaluateTaylor(
                   coeff, BoundaryClosure::taylorOrder(accuracyOrder),
                   BoundaryClosure::ghostDistance(layer, h)));
        return;
    }

    /// 实边界点：直接赋 BC 值，并向该法向的虚胞层广播。
    setter(i, j, k, bcValue);

    auto broadcast = [&](int gi, int gj, int gk,
                         int ri, int rj, int rk) {
        /// 确保参考内点所有坐标都在物理域内。
        int refI = ri, refJ = rj, refK = rk;
        nearestInteriorCell(field, ri, rj, rk, refI, refJ, refK);

        std::array<double, BoundaryClosure::kMaxAccuracyOrder> coeff{};
        BoundaryClosure::buildMultiDimScalarTaylorCoefficients(
            field, refI, refJ, refK, normal, accuracyOrder,
            getter, bcName, coeff);
        coeff[0] = bcValue;

        const int layer =
            BoundaryClosure::ghostLayer(gi, gj, gk, refI, refJ, refK, normal);
        if (layer < 1 || layer > field.NG()) {
            BoundaryClosure::fatalClosure(bcName, gi, gj, gk,
                                  "invalid ghost layer for fixed-value closure.");
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

void applyScalar(Field& field, int i, int j, int k,
                 double bcValue, int vIdx, int accuracyOrder) {
    auto setter = [&](int gi, int gj, int gk, double value) {
        field(gi, gj, gk, vIdx) = value;
    };
    writeFixedValueGhosts(field, i, j, k, bcValue, accuracyOrder,
                          "fixedValue",
                          setter, BoundaryClosure::FieldVariableGetter{vIdx});
}

void applyVector3(Field& field, int i, int j, int k,
                  const SF::Vector3& bcValue, int vIdx,
                  int accuracyOrder) {
    if (vIdx != RU) {
        throw std::runtime_error(
            "ILW fixedValue Vector3 boundary is only implemented for velocity/RU.");
    }
    const double fixed[3] = {bcValue.x, bcValue.y, bcValue.z};
    for (int comp = 0; comp < 3; ++comp) {
        auto setter = [&](int gi, int gj, int gk, double value) {
            const double rho =
                requireDensity(field, gi, gj, gk,
                               "ILW fixedValue velocity write");
            field(gi, gj, gk, vIdx + comp) = rho * value;
        };
        writeFixedValueGhosts(field, i, j, k, fixed[comp], accuracyOrder,
                              "fixedValueVector3",
                              setter,
                              PrimitiveVelocityComponentGetter{comp});
    }
}

} // namespace FixedValue
} // namespace ILW
} // namespace Boundary
} // namespace SF
