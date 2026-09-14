/// @file SF_zeroGradient.cpp
/// @brief 线性 ghost 层的 zeroGradient 标量/向量实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_zeroGradient.h"
#include "core/mesh/SF_dimension.h"

#include <algorithm>

namespace SF {
namespace Boundary {
namespace Algebraic {
namespace ZeroGradient {
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

void nearestInsideSample(const Field& field,
                         int i, int j, int k, int axis,
                         int& si, int& sj, int& sk) {
    const int ng = field.NG();
    si = i;
    sj = j;
    sk = k;
    if (axis == 0) {
        if (i == ng && field.NX() > 1) si = i + 1;
        else if (i == field.NX() + ng - 1 && field.NX() > 1) si = i - 1;
    } else if (axis == 1) {
        if (j == ng && field.NY() > 1) sj = j + 1;
        else if (j == field.NY() + ng - 1 && field.NY() > 1) sj = j - 1;
    } else if (axis == 2) {
        if (k == ng && field.NZ() > 2) sk = k + 1;
        else if (k == field.NZ() + ng - 1 && field.NZ() > 2) sk = k - 1;
    }
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

} // namespace

void applyScalar(Field& field, int i, int j, int k, int axis, int vIdx) {
    if (isGhostCell(field, i, j, k)) {
        int ri = 0, rj = 0, rk = 0;
        nearestInteriorCell(field, i, j, k, ri, rj, rk);
        field(i, j, k, vIdx) = field(ri, rj, rk, vIdx);
        return;
    }

    int si = i, sj = j, sk = k;
    nearestInsideSample(field, i, j, k, axis, si, sj, sk);
    field(i, j, k, vIdx) = field(si, sj, sk, vIdx);

    auto broadcast = [&](int gi, int gj, int gk,
                         int ri, int rj, int rk) {
        field(gi, gj, gk, vIdx) = field(ri, rj, rk, vIdx);
    };
    forBoundaryGhostsAlongAxis(field, i, j, k, axis, broadcast);
}

void applyVector3(Field& field, int i, int j, int k, int axis, int vIdx) {
    if (isGhostCell(field, i, j, k)) {
        int ri = 0, rj = 0, rk = 0;
        nearestInteriorCell(field, i, j, k, ri, rj, rk);
        field.setVal<SF::Vector3>(
            i, j, k, vIdx,
            field.getVal<SF::Vector3>(ri, rj, rk, vIdx));
        return;
    }

    int si = i, sj = j, sk = k;
    nearestInsideSample(field, i, j, k, axis, si, sj, sk);
    field.setVal<SF::Vector3>(
        i, j, k, vIdx,
        field.getVal<SF::Vector3>(si, sj, sk, vIdx));

    auto broadcast = [&](int gi, int gj, int gk,
                         int ri, int rj, int rk) {
        field.setVal<SF::Vector3>(
            gi, gj, gk, vIdx,
            field.getVal<SF::Vector3>(ri, rj, rk, vIdx));
    };
    forBoundaryGhostsAlongAxis(field, i, j, k, axis, broadcast);
}

} // namespace ZeroGradient
} // namespace Algebraic
} // namespace Boundary
} // namespace SF
