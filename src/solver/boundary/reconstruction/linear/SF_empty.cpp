/// @file SF_empty.cpp
/// @brief 线性重建中的 empty 非活动维边界实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_empty.h"

namespace SF {
namespace Boundary {
namespace Algebraic {
namespace Empty {
namespace {

void offsetForAxis(int axis, int& di, int& dj, int& dk) {
    di = dj = dk = 0;
    if (axis == 0) di = 1;
    else if (axis == 1) dj = 1;
    else dk = 1;
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

} // namespace

void applyScalar(Field& field, int i, int j, int k, int vIdx, int axis) {
    int di = 0, dj = 0, dk = 0;
    offsetForAxis(axis, di, dj, dk);

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
    int di = 0, dj = 0, dk = 0;
    offsetForAxis(axis, di, dj, dk);

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
} // namespace Algebraic
} // namespace Boundary
} // namespace SF
