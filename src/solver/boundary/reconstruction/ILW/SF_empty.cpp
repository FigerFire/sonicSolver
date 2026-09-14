/// @file SF_empty.cpp
/// @brief ILW 管线中的 empty 非活动维边界处理。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_empty.h"

#include "SF_ILW.h"
#include "core/mesh/SF_dimension.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>

namespace SF {
namespace Boundary {
namespace ILW {
namespace Empty {
namespace {

void fatalInvalidEmptyAxis(int axis) {
    std::cerr << "[SF FATAL] ILW empty boundary received invalid axis="
              << axis << ". EMPTY zones must resolve to one structured "
              << "boundary direction." << std::endl;
    std::exit(1);
}

void offsetForAxis(int axis, int& di, int& dj, int& dk) {
    di = dj = dk = 0;
    if (axis == 0) di = 1;
    else if (axis == 1) dj = 1;
    else if (axis == 2) dk = 1;
    else fatalInvalidEmptyAxis(axis);
}

bool isLowerBoundary(const Field& field, int i, int j, int k, int axis) {
    const int ng = field.NG();
    return (axis == 0 && i == ng) ||
           (axis == 1 && j == ng) ||
           (axis == 2 && k == ng);
}

bool isUpperBoundary(const Field& field, int i, int j, int k, int axis) {
    const int ng = field.NG();
    return (axis == 0 && i == field.NX() + ng - 1) ||
           (axis == 1 && j == field.NY() + ng - 1) ||
           (axis == 2 && k == field.NZ() + ng - 1);
}

int clampToInteriorAxis(const Field& field, int coord, int axis) {
    const int ng = field.NG();
    const int n = (axis == 0) ? field.NX()
                : (axis == 1) ? field.NY()
                              : field.NZ();
    return std::max(ng, std::min(coord, ng + n - 1));
}

void nearestEmptyInteriorCell(const Field& field,
                              int i, int j, int k,
                              int& ri, int& rj, int& rk) {
    ri = clampToInteriorAxis(field, i, 0);
    rj = clampToInteriorAxis(field, j, 1);
    rk = clampToInteriorAxis(field, k, 2);
}

void ensureEmptyDimension(int axis) {
    if (axis < 0 || axis > 2) fatalInvalidEmptyAxis(axis);

    const int inactive = Math::singleInactiveDirectionIndex();
    if (inactive >= 0 && inactive != axis) {
        std::cerr << "[SF FATAL] ILW empty boundary requests inactive axis "
                  << axis << ", but current dimension state already closed axis "
                  << inactive << ". Multiple EMPTY directions are not supported."
                  << std::endl;
        std::exit(1);
    }

    if (Math::activeDimensionCount() == 3) {
        Math::deactivateDirection(axis);
        static bool emitted[3] = {false, false, false};
        if (!emitted[axis]) {
            std::cerr << "[SF ILW] EMPTY boundary deactivates direction "
                      << axis << "; subsequent ILW boundary closure uses "
                      << Math::activeDimensionCount() << "D active space."
                      << std::endl;
            emitted[axis] = true;
        }
    }
}

} // namespace

void applyScalar(Field& field, int i, int j, int k, int vIdx, int axis) {
    ensureEmptyDimension(axis);

    int di = 0, dj = 0, dk = 0;
    offsetForAxis(axis, di, dj, dk);

    if (ILW::isGhostCell(field, i, j, k)) {
        int ri = 0, rj = 0, rk = 0;
        nearestEmptyInteriorCell(field, i, j, k, ri, rj, rk);
        field(i, j, k, vIdx) = field(ri, rj, rk, vIdx);
        return;
    }

    if (isLowerBoundary(field, i, j, k, axis)) {
        for (int b = 1; b <= field.NG(); ++b) {
            field(i - b * di, j - b * dj, k - b * dk, vIdx) =
                field(i, j, k, vIdx);
        }
    }
    if (isUpperBoundary(field, i, j, k, axis)) {
        for (int b = 1; b <= field.NG(); ++b) {
            field(i + b * di, j + b * dj, k + b * dk, vIdx) =
                field(i, j, k, vIdx);
        }
    }
}

void applyVector3(Field& field, int i, int j, int k, int vIdx, int axis) {
    ensureEmptyDimension(axis);

    int di = 0, dj = 0, dk = 0;
    offsetForAxis(axis, di, dj, dk);

    if (ILW::isGhostCell(field, i, j, k)) {
        int ri = 0, rj = 0, rk = 0;
        nearestEmptyInteriorCell(field, i, j, k, ri, rj, rk);
        field.setVal<SF::Vector3>(
            i, j, k, vIdx,
            field.getVal<SF::Vector3>(ri, rj, rk, vIdx));
        return;
    }

    if (isLowerBoundary(field, i, j, k, axis)) {
        for (int b = 1; b <= field.NG(); ++b) {
            field.setVal<SF::Vector3>(
                i - b * di, j - b * dj, k - b * dk, vIdx,
                field.getVal<SF::Vector3>(i, j, k, vIdx));
        }
    }
    if (isUpperBoundary(field, i, j, k, axis)) {
        for (int b = 1; b <= field.NG(); ++b) {
            field.setVal<SF::Vector3>(
                i + b * di, j + b * dj, k + b * dk, vIdx,
                field.getVal<SF::Vector3>(i, j, k, vIdx));
        }
    }
}

} // namespace Empty
} // namespace ILW
} // namespace Boundary
} // namespace SF
