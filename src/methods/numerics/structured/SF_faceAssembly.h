#pragma once

/// @file SF_faceAssembly.h
/// @brief 数值面通量存储及结构方向残差装配。

#include "methods/numerics/structured/SF_faceGeometry.h"
#include "core/flux/SF_flux.h"
#include "core/residual/SF_residual.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace SF::Math {

/// @brief GlobalFace 仅由 canonical owner 执行数值通量核。
inline bool shouldCalculateFaceFlux(
        const Field& field, Dir d, int i, int j, int k) {
    return !field.hasCanonicalFaceMetrics(
               static_cast<int>(d), i, j, k)
        || field.isCanonicalFaceOwner(
               static_cast<int>(d), i, j, k);
}

inline void storeFlux(FluxField& fluxField, const Field& field,
                      int i, int j, int k, Dir d,
                      const double flux[5]) {
    for (int v = 0; v < 5; ++v) {
        if (!std::isfinite(flux[v])) {
            std::cerr << "[SF FATAL] non-finite conservative face flux at face ("
                      << i << "," << j << "," << k << "), dir=" << d
                      << ", component=" << v << ", value=" << flux[v]
                      << std::endl;
            std::exit(1);
        }
        fluxField(static_cast<std::size_t>(d) * field.TotalSize()
                      + field.getIdx(i, j, k), v) = flux[v];
    }
}

inline void storeFlux(FluxField& fluxField, const Field& field,
                      int i, int j, int k, Dir d,
                      const std::vector<double>& flux) {
    if (static_cast<int>(flux.size()) != field.NVar()) {
        std::cerr << "[SF FATAL] FluidStateModel face flux size mismatch."
                  << std::endl;
        std::exit(1);
    }
    for (int v = 0; v < field.NVar(); ++v) {
        if (!std::isfinite(flux[static_cast<size_t>(v)])) {
            std::cerr << "[SF FATAL] non-finite FluidStateModel face flux."
                      << std::endl;
            std::exit(1);
        }
        fluxField(static_cast<std::size_t>(d) * field.TotalSize()
                      + field.getIdx(i, j, k), v) =
            flux[static_cast<size_t>(v)];
    }
}

inline void assembleConvectiveFluxResidual(
        Field& field, const FluxField& fluxField, Residual& residual) {
    auto assembleDirection = [&](Dir d) {
        forFaces(field, d, [&](int i, int j, int k) {
            for (int v = 0; v < field.NVar(); ++v) {
                const double flux =
                    fluxField(static_cast<std::size_t>(d) * field.TotalSize()
                              + field.getIdx(i, j, k), v);
                if (!std::isfinite(flux)) {
                    std::cerr
                        << "[SF FATAL] non-finite stored face flux while "
                        << "assembling residual at face (" << i << ","
                        << j << "," << k << "), dir=" << d
                        << ", component=" << v << std::endl;
                    std::exit(1);
                }
                if (d == XI) residual.x(i, j, k, v) = flux;
                else if (d == ETA) residual.y(i, j, k, v) = flux;
                else residual.z(i, j, k, v) = flux;
            }
        });
    };

    assembleDirection(XI);
    assembleDirection(ETA);
    assembleDirection(ZETA);
}

} // namespace SF::Math
