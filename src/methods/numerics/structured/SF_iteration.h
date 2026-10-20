#pragma once

/// @file SF_iteration.h
/// @brief 结构网格方向和 cell/face 遍历工具。

#include "core/mesh/SF_dimension.h"
#include "core/field/SF_field.h"

#include <cmath>

namespace SF::Math {

enum Dir { XI = 0, ETA = 1, ZETA = 2 };

inline bool isDirectionActive(Dir d) {
    return isDirectionActiveIndex(static_cast<int>(d));
}

inline void dirOffset(Dir d, int& di, int& dj, int& dk) {
    di = dj = dk = 0;
    if (d == XI) di = 1;
    else if (d == ETA) dj = 1;
    else dk = 1;
}

template <typename Func>
inline void forInterior(const Field& field, Func&& body) {
    const int ng = field.NG();
    const int nx = field.NX(), ny = field.NY(), nz = field.NZ();
    for (int k = ng; k < nz + ng; ++k)
        for (int j = ng; j < ny + ng; ++j)
            for (int i = ng; i < nx + ng; ++i)
                body(i, j, k);
}

template <typename Func>
inline void forSolvedCells(const Field& field, Func&& body) {
    const int ng = field.NG();
    const int nx = field.NX(), ny = field.NY(), nz = field.NZ();
    for (int k = ng; k < nz + ng; ++k)
        for (int j = ng; j < ny + ng; ++j)
            for (int i = ng; i < nx + ng; ++i) {
                if (field.isSolverBoundaryPoint(i, j, k)) continue;
                body(i, j, k);
            }
}

inline bool isSolvedFluidCell(const Field& field, int i, int j, int k) {
    return field.CellFlag(i, j, k) == FLUID_CELL;
}

template <typename Func>
inline void forFluidInterior(const Field& field, Func&& body) {
    forSolvedCells(field, [&](int i, int j, int k) {
        if (!isSolvedFluidCell(field, i, j, k)) return;
        body(i, j, k);
    });
}

template <typename Func>
inline void forFaces(const Field& field, Dir d, Func&& body) {
    if (!isDirectionActive(d)) return;

    const int ng = field.NG();
    const int nx = field.NX(), ny = field.NY(), nz = field.NZ();
    if (d == XI) {
        for (int k = ng; k < nz + ng; ++k)
            for (int j = ng; j < ny + ng; ++j)
                for (int i = ng - 1; i < nx + ng; ++i)
                    body(i, j, k);
    } else if (d == ETA) {
        for (int k = ng; k < nz + ng; ++k)
            for (int j = ng - 1; j < ny + ng; ++j)
                for (int i = ng; i < nx + ng; ++i)
                    body(i, j, k);
    } else {
        for (int k = ng - 1; k < nz + ng; ++k)
            for (int j = ng; j < ny + ng; ++j)
                for (int i = ng; i < nx + ng; ++i)
                    body(i, j, k);
    }
}

inline double vecMag(double x, double y, double z) {
    return std::sqrt(x * x + y * y + z * z);
}

} // namespace SF::Math
