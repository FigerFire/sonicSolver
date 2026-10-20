#pragma once

/// @file SF_residualAssembly.h
/// @brief 已装配结构方向通量的空间残差和显式状态更新。

#include "methods/numerics/structured/SF_iteration.h"
#include "methods/numerics/SF_utility.h"
#include "core/residual/SF_residual.h"

namespace SF::Math {

inline double localSpatialResidual(const Field& field, const Residual& residual,
                                   int i, int j, int k, int variable) {
    double transformedFluxDivergence = 0.0;
    if (isDirectionActive(XI)) {
        transformedFluxDivergence +=
            residual.x(i, j, k, variable)
            - residual.x(i - 1, j, k, variable);
    }
    if (isDirectionActive(ETA)) {
        transformedFluxDivergence +=
            residual.y(i, j, k, variable)
            - residual.y(i, j - 1, k, variable);
    }
    if (isDirectionActive(ZETA)) {
        transformedFluxDivergence +=
            residual.z(i, j, k, variable)
            - residual.z(i, j, k - 1, variable);
    }
    return field.Jac(i, j, k) * transformedFluxDivergence
         - residual.source(i, j, k, variable);
}

/// @brief 返回 Runtime 汇总后的 GlobalDof 残差；普通点返回本地强残差。
inline double spatialResidual(const Field& field, const Residual& residual,
                              int i, int j, int k, int variable) {
    if (residual.hasGlobal(i, j, k)) {
        return residual.global(i, j, k, variable);
    }
    return localSpatialResidual(field, residual, i, j, k, variable);
}

/// @brief 在 MPI 无关的数值层冻结本地强残差，供 Runtime 后续积分装配。
inline void stageLocalSpatialResidual(Field& field, Residual& residual) {
    residual.clearGlobal();
    forFluidInterior(field, [&](int i, int j, int k) {
        for (int v = 0; v < field.NVar(); ++v) {
            residual.stageLocal(
                i, j, k, v,
                localSpatialResidual(field, residual, i, j, k, v));
        }
    });
}

inline void applyDivergence(Field& field, const Residual& residual, double dt) {
    forFluidInterior(field, [&](int i, int j, int k) {
        for (int v = 0; v < field.NVar(); ++v) {
            field(i, j, k, v) -= dt * spatialResidual(field, residual, i, j, k, v);
        }
        if (field.NVar() == 5) {
            Numerics::requirePhysicalState(
                "post-update conservative state", field, i, j, k);
        }
    });
}

} // namespace SF::Math
