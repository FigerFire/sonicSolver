/// @file SF_symmetry.cpp
/// @brief 线性对称边界的法向速度反射实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_symmetry.h"

#include "SF_zeroGradient.h"

#include <algorithm>

namespace SF {
namespace Boundary {
namespace Algebraic {
namespace Symmetry {
namespace {

void nearestInteriorCell(const Field& field,
                         int i, int j, int k,
                         int& ri, int& rj, int& rk) {
    const int ng = field.NG();
    ri = std::max(ng, std::min(i, field.NX() + ng - 1));
    rj = std::max(ng, std::min(j, field.NY() + ng - 1));
    rk = std::max(ng, std::min(k, field.NZ() + ng - 1));
}

void reflectNormal(SF::Vector3& value, int axis) {
    if (axis == 0) value.x = -value.x;
    else if (axis == 1) value.y = -value.y;
    else if (axis == 2) value.z = -value.z;
}

void zeroNormal(SF::Vector3& value, int axis) {
    if (axis == 0) value.x = 0.0;
    else if (axis == 1) value.y = 0.0;
    else if (axis == 2) value.z = 0.0;
}

template <typename Fn>
void forBoundaryGhostsAlongAxis(Field& field,
                                int i, int j, int k,
                                int axis,
                                Fn&& fn) {
    const int ng = field.NG();
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

} // namespace

void applyScalar(Field& field, int i, int j, int k, int axis, int vIdx) {
    ZeroGradient::applyScalar(field, i, j, k, axis, vIdx);
}

void applyVector3(Field& field, int i, int j, int k, int axis, int vIdx) {
    int ri = 0, rj = 0, rk = 0;
    nearestInteriorCell(field, i, j, k, ri, rj, rk);

    SF::Vector3 value =
        field.getVal<SF::Vector3>(ri, rj, rk, vIdx);
    const int ng = field.NG();
    if (i < ng || i >= field.NX() + ng ||
        j < ng || j >= field.NY() + ng ||
        k < ng || k >= field.NZ() + ng) {
        reflectNormal(value, axis);
        field.setVal<SF::Vector3>(i, j, k, vIdx, value);
        return;
    }

    zeroNormal(value, axis);
    field.setVal<SF::Vector3>(i, j, k, vIdx, value);

    auto reflectGhost = [&](int gi, int gj, int gk,
                            int bi, int bj, int bk) {
        SF::Vector3 ghostValue =
            field.getVal<SF::Vector3>(bi, bj, bk, vIdx);
        reflectNormal(ghostValue, axis);
        field.setVal<SF::Vector3>(gi, gj, gk, vIdx, ghostValue);
    };
    forBoundaryGhostsAlongAxis(field, i, j, k, axis, reflectGhost);
}

} // namespace Symmetry
} // namespace Algebraic
} // namespace Boundary
} // namespace SF
