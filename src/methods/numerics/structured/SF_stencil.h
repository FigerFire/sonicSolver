#pragma once

/// @file SF_stencil.h
/// @brief 结构方向 stencil 常量、深度检查和状态提取。

#include "methods/numerics/structured/SF_iteration.h"

#include <iostream>

namespace SF::Math {

inline constexpr int WENO3_OFFSETS[4] = {-1, 0, 1, 2};
inline constexpr int WENO3_STENCIL = 4;
inline constexpr int WENO5_OFFSETS[6] = {-2, -1, 0, 1, 2, 3};
inline constexpr int WENO5_STENCIL = 6;
inline constexpr int WENO7_OFFSETS[8] = {-3, -2, -1, 0, 1, 2, 3, 4};
inline constexpr int WENO7_STENCIL = 8;

inline constexpr int minGhostWENO3 = 2;
inline constexpr int minGhostWENO5 = 3;
inline constexpr int minGhostWENO7 = 4;

inline bool checkGhostDepth(const Field& field,
                            int requiredGhost,
                            const char* schemeName) {
    const int ng = field.NG();
    if (ng < requiredGhost) {
        std::cerr << "[SF FATAL] " << schemeName
                  << " requires ghost >= " << requiredGhost
                  << ", but field.NG() = " << ng << std::endl;
        return false;
    }
    return true;
}

inline void extractStencil(const Field& field,
                           int i, int j, int k, Dir d,
                           const int* offsets, int nStencil,
                           double* output) {
    int di = 0, dj = 0, dk = 0;
    dirOffset(d, di, dj, dk);
    for (int s = 0; s < nStencil; ++s) {
        const int ii = i + offsets[s] * di;
        const int jj = j + offsets[s] * dj;
        const int kk = k + offsets[s] * dk;
        double* row = output + s * 5;
        for (int v = 0; v < 5; ++v) row[v] = field(ii, jj, kk, v);
    }
}

} // namespace SF::Math
