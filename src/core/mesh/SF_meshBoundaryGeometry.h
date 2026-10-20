/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_meshBoundaryGeometry.h
/// @brief 中立结构网格边界拓扑和几何查询。

#include "core/field/SF_field.h"
#include "core/mesh/SF_dimension.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace SF::StructuredMesh::BoundaryGeometry {

/// @brief 边界 set 的结构化拓扑信息。
struct BoundarySetInfo {
    bool valid = false;
    int axis = -1;
    std::array<double, 3> normal{0.0, 0.0, 1.0};
};

inline int axisSize(const Field& field, int axis) {
    if (axis == 0) return field.NX();
    if (axis == 1) return field.NY();
    return field.NZ();
}

inline int coordinateIndexForAxis(int i, int j, int k, int axis) {
    if (axis == 0) return i;
    if (axis == 1) return j;
    return k;
}

inline bool onAxisBoundary(const Field& field,
                           int i, int j, int k,
                           int axis) {
    const int lo = field.NG();
    const int hi = lo + axisSize(field, axis) - 1;
    const int coord = coordinateIndexForAxis(i, j, k, axis);
    return coord == lo || coord == hi;
}

inline bool isPhysicalPoint(const Field& field, int i, int j, int k) {
    const int ng = field.NG();
    return i >= ng && i < ng + field.NX()
        && j >= ng && j < ng + field.NY()
        && k >= ng && k < ng + field.NZ();
}

inline std::array<double, 3> normalized(std::array<double, 3> v) {
    const double mag = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (mag <= 1.0e-14) return {0.0, 0.0, 1.0};
    return {v[0] / mag, v[1] / mag, v[2] / mag};
}

inline double dot(const std::array<double, 3>& a,
                  const std::array<double, 3>& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline std::array<double, 3> cross(const std::array<double, 3>& a,
                                   const std::array<double, 3>& b) {
    return {
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    };
}

inline std::array<double, 3> point(const Field& field, int i, int j, int k) {
    return {field.X(i, j, k), field.Y(i, j, k), field.Z(i, j, k)};
}

inline std::array<double, 3> subtract(const std::array<double, 3>& a,
                                      const std::array<double, 3>& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

inline std::array<double, 3> add(const std::array<double, 3>& a,
                                 const std::array<double, 3>& b) {
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

inline std::array<double, 3> estimateAxisNormal(const Field& field, int axis) {
    const int ng = field.NG();
    const int i0 = ng;
    const int j0 = ng;
    const int k0 = ng;
    const int i1 = ng + field.NX() - 1;
    const int j1 = ng + field.NY() - 1;
    const int k1 = ng + field.NZ() - 1;

    if (axisSize(field, axis) >= 2) {
        std::array<double, 3> sum{0.0, 0.0, 0.0};
        int n = 0;
        if (axis == 0) {
            for (int k = k0; k <= k1; ++k)
                for (int j = j0; j <= j1; ++j) {
                    sum = add(sum, subtract(point(field, i0 + 1, j, k),
                                            point(field, i0, j, k)));
                    ++n;
                }
        } else if (axis == 1) {
            for (int k = k0; k <= k1; ++k)
                for (int i = i0; i <= i1; ++i) {
                    sum = add(sum, subtract(point(field, i, j0 + 1, k),
                                            point(field, i, j0, k)));
                    ++n;
                }
        } else {
            for (int j = j0; j <= j1; ++j)
                for (int i = i0; i <= i1; ++i) {
                    sum = add(sum, subtract(point(field, i, j, k0 + 1),
                                            point(field, i, j, k0)));
                    ++n;
                }
        }
        if (n > 0) return normalized(sum);
    }

    const int a0 = (axis + 1) % 3;
    const int a1 = (axis + 2) % 3;
    auto shifted = [&](int i, int j, int k, int a) {
        if (a == 0 && field.NX() > 1) i = std::min(i + 1, i1);
        if (a == 1 && field.NY() > 1) j = std::min(j + 1, j1);
        if (a == 2 && field.NZ() > 1) k = std::min(k + 1, k1);
        return point(field, i, j, k);
    };

    const auto p0 = point(field, i0, j0, k0);
    const auto ta = subtract(shifted(i0, j0, k0, a0), p0);
    const auto tb = subtract(shifted(i0, j0, k0, a1), p0);
    return normalized(cross(ta, tb));
}

inline BoundarySetInfo analyzeBoundarySet(const Field& field,
                                          const std::string& name) {
    const auto& sets = field.getAllSets();
    auto it = sets.find(name);
    if (it == sets.end() || it->second.empty()) return {};

    std::vector<int> candidates;
    for (int axis = 0; axis < 3; ++axis) {
        bool allOnBoundary = true;
        int nPhysical = 0;
        for (int idx : it->second) {
            int i = 0, j = 0, k = 0;
            field.getIJK(idx, i, j, k);
            if (!isPhysicalPoint(field, i, j, k)) continue;
            ++nPhysical;
            if (!onAxisBoundary(field, i, j, k, axis)) {
                allOnBoundary = false;
                break;
            }
        }
        if (nPhysical > 0 && allOnBoundary) candidates.push_back(axis);
    }

    if (candidates.size() != 1) return {};

    BoundarySetInfo info;
    info.valid = true;
    info.axis = candidates.front();
    info.normal = estimateAxisNormal(field, info.axis);
    return info;
}

inline int boundaryAxisForSet(const Field& field, const std::string& name) {
    const BoundarySetInfo info = analyzeBoundarySet(field, name);
    if (!info.valid) {
        std::cerr << "[SF FATAL] Boundary set '" << name
                  << "' is not located on one structured boundary direction."
                  << std::endl;
        std::exit(1);
    }
    return info.axis;
}

inline int activeBoundaryAxisForSet(const Field& field,
                                    const std::string& name) {
    const auto& sets = field.getAllSets();
    auto it = sets.find(name);
    if (it == sets.end() || it->second.empty()) {
        std::cerr << "[SF FATAL] Boundary set '" << name
                  << "' is not a mesh set." << std::endl;
        std::exit(1);
    }

    std::vector<int> candidates;
    for (int axis = 0; axis < 3; ++axis) {
        if (!Math::isDirectionActiveIndex(axis)) continue;
        bool allOnBoundary = true;
        int nPhysical = 0;
        for (int idx : it->second) {
            int i = 0, j = 0, k = 0;
            field.getIJK(idx, i, j, k);
            if (!isPhysicalPoint(field, i, j, k)) continue;
            ++nPhysical;
            if (!onAxisBoundary(field, i, j, k, axis)) {
                allOnBoundary = false;
                break;
            }
        }
        if (nPhysical > 0 && allOnBoundary) candidates.push_back(axis);
    }

    if (candidates.size() != 1) {
        std::cerr << "[SF FATAL] Boundary set '" << name
                  << "' is not located on one active structured boundary "
                  << "direction." << std::endl;
        std::exit(1);
    }
    return candidates.front();
}

} // namespace SF::StructuredMesh::BoundaryGeometry
