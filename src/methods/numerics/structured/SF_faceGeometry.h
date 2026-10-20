#pragma once

/// @file SF_faceGeometry.h
/// @brief 结构网格守恒面度量访问。

#include "methods/numerics/structured/SF_canonicalFace.h"
#include "methods/numerics/structured/SF_iteration.h"

namespace SF::Math {

inline double metricCofactor(double metric, double inverseJacobian) {
    return Numerics::CanonicalFace::metricCofactor(metric, inverseJacobian);
}

inline void faceMetrics(const Field& field,
                        int i, int j, int k, Dir d,
                        double metrics[4]) {
    const auto geometry = Numerics::CanonicalFace::geometry(
        field, static_cast<int>(d), i, j, k);
    metrics[0] = geometry.cofactor[0];
    metrics[1] = geometry.cofactor[1];
    metrics[2] = geometry.cofactor[2];
    metrics[3] = geometry.inverseJacobian;
}

} // namespace SF::Math
