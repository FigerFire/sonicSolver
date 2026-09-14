/// @file SF_symmetry.cpp
/// @brief ILW 对称/滑移壁面的法向特征闭合。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_symmetry.h"

#include "SF_ILW.h"
#include "core/mesh/SF_dimension.h"
#include "SF_zeroGradient.h"

namespace SF {
namespace Boundary {
namespace ILW {
namespace Symmetry {
namespace {

int activeBoundaryAxisForPoint(const Field& field, int i, int j, int k) {
    const int ng = field.NG();
    if (Math::isDirectionActiveIndex(0) &&
        (i == ng || i == field.NX() + ng - 1 ||
         i < ng || i >= field.NX() + ng)) {
        return 0;
    }
    if (Math::isDirectionActiveIndex(1) &&
        (j == ng || j == field.NY() + ng - 1 ||
         j < ng || j >= field.NY() + ng)) {
        return 1;
    }
    if (Math::isDirectionActiveIndex(2) &&
        (k == ng || k == field.NZ() + ng - 1 ||
         k < ng || k >= field.NZ() + ng)) {
        return 2;
    }
    return -1;
}

void zeroNormalComponent(SF::Vector3& value, int axis) {
    if (axis == 0) value.x = 0.0;
    else if (axis == 1) value.y = 0.0;
    else if (axis == 2) value.z = 0.0;
}

void reflectNormalComponentByAxis(SF::Vector3& value, int axis) {
    if (axis == 0) value.x = -value.x;
    else if (axis == 1) value.y = -value.y;
    else if (axis == 2) value.z = -value.z;
}

} // namespace

void applyScalar(Field& field, int i, int j, int k, int vIdx,
                 int accuracyOrder) {
    ZeroGradient::applyScalar(field, i, j, k, vIdx, accuracyOrder);
}

void applyVector3(Field& field, int i, int j, int k, int vIdx,
                  int accuracyOrder) {
    if (isGhostCell(field, i, j, k)) {
        ZeroGradient::applyVector3(field, i, j, k, vIdx, accuracyOrder);
        SF::Vector3 value = field.getVal<SF::Vector3>(i, j, k, vIdx);
        const int axis = activeBoundaryAxisForPoint(field, i, j, k);
        reflectNormalComponentByAxis(value, axis);
        field.setVal<SF::Vector3>(i, j, k, vIdx, value);
        return;
    }

    ZeroGradient::applyVector3(field, i, j, k, vIdx, accuracyOrder);

    const int boundaryAxis = activeBoundaryAxisForPoint(field, i, j, k);
    if (boundaryAxis >= 0) {
        SF::Vector3 boundaryValue =
            field.getVal<SF::Vector3>(i, j, k, vIdx);
        zeroNormalComponent(boundaryValue, boundaryAxis);
        field.setVal<SF::Vector3>(i, j, k, vIdx, boundaryValue);
    }

    auto reflectGhost = [&](int gi, int gj, int gk,
                            int ri, int rj, int rk) {
        (void)ri;
        (void)rj;
        (void)rk;
        SF::Vector3 value =
            field.getVal<SF::Vector3>(gi, gj, gk, vIdx);
        reflectNormalComponentByAxis(value, boundaryAxis);
        field.setVal<SF::Vector3>(gi, gj, gk, vIdx, value);
    };
    forBoundaryGhostsAlongAxis(field, i, j, k, boundaryAxis, reflectGhost);
}

} // namespace Symmetry
} // namespace ILW
} // namespace Boundary
} // namespace SF
