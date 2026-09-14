/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_boundary.h
/// @brief 物理边界条件统一调度入口。

#include "reconstruction/SF_reconstruction.h"
#include "reconstruction/ILW/SF_boundaryClosure.h"
#include "SF_boundaryGeometry.h"
#include "core/mesh/SF_dimension.h"
#include "SF_thermodynamicClosure.h"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <type_traits>

namespace SF {
namespace Boundary {

template <typename T>
struct AlwaysFalse : std::false_type {};

/// @brief 拷贝普通虚胞的网格度规。
void copyMetrics(Field& field, int i, int j, int k,
                 int ii, int jj, int kk);

/// @brief 计算普通物理边界虚胞的镜像实胞索引。
void getMirrorIJK(int i, int j, int k,
                  int ng, int nx, int ny, int nz,
                  int& ii, int& jj, int& kk);

/// @brief 判断索引是否位于普通物理边界虚胞层。
inline bool isBoundaryGhostCell(const Field& field, int i, int j, int k) {
    const int ng = field.NG();
    return (i < ng || i >= field.NX() + ng ||
            j < ng || j >= field.NY() + ng ||
            k < ng || k >= field.NZ() + ng);
}

/// @brief 将普通边界索引夹到最近的物理域实胞。
inline void nearestBoundaryInteriorCell(const Field& field,
                                        int i, int j, int k,
                                        int& ri, int& rj, int& rk) {
    const int ng = field.NG();
    ri = std::max(ng, std::min(i, field.NX() + ng - 1));
    rj = std::max(ng, std::min(j, field.NY() + ng - 1));
    rk = std::max(ng, std::min(k, field.NZ() + ng - 1));
}

/// @brief 对实边界点外侧虚胞层执行回调。
///
/// 只填充活跃方向的外层虚胞；被 EMPTY 关闭的方向不在此处填充，
/// 避免与 EMPTY 边界条件产生角点冲突。
template <typename Fn>
void forBoundaryGhosts(Field& field, int i, int j, int k, Fn&& fn) {
    const int ng = field.NG();
    if (Math::isDirectionActiveIndex(0)) {
        if (i == ng) {
            for (int b = 0; b < ng; ++b) fn(b, j, k, i, j, k);
        }
        if (i == field.NX() + ng - 1) {
            for (int b = 1; b <= ng; ++b) fn(i + b, j, k, i, j, k);
        }
    }
    if (Math::isDirectionActiveIndex(1)) {
        if (j == ng) {
            for (int b = 0; b < ng; ++b) fn(i, b, k, i, j, k);
        }
        if (j == field.NY() + ng - 1) {
            for (int b = 1; b <= ng; ++b) fn(i, j + b, k, i, j, k);
        }
    }
    if (Math::isDirectionActiveIndex(2)) {
        if (k == ng) {
            for (int b = 0; b < ng; ++b) fn(i, j, b, i, j, k);
        }
        if (k == field.NZ() + ng - 1) {
            for (int b = 1; b <= ng; ++b) fn(i, j, k + b, i, j, k);
        }
    }
}

/// @brief 只沿指定普通物理边界法向填充外侧虚胞层。
///
/// 该函数用于按 patch 法向应用边界条件，避免角点处一个 patch 的 BC
/// 写入另一个 patch 的 ghost 层。
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

/// @brief 遍历某个计算方向两端的物理面点。
template <typename Fn>
void forAxisBoundaryPoints(Field& field, int axis, Fn&& fn) {
    const int ng = field.NG();
    const int i0 = ng;
    const int i1 = field.NX() + ng - 1;
    const int j0 = ng;
    const int j1 = field.NY() + ng - 1;
    const int k0 = ng;
    const int k1 = field.NZ() + ng - 1;

    if (axis == 0) {
        for (int k = k0; k <= k1; ++k)
            for (int j = j0; j <= j1; ++j) {
                fn(i0, j, k);
                if (i1 != i0) fn(i1, j, k);
            }
    } else if (axis == 1) {
        for (int k = k0; k <= k1; ++k)
            for (int i = i0; i <= i1; ++i) {
                fn(i, j0, k);
                if (j1 != j0) fn(i, j1, k);
            }
    } else if (axis == 2) {
        for (int j = j0; j <= j1; ++j)
            for (int i = i0; i <= i1; ++i) {
                fn(i, j, k0);
                if (k1 != k0) fn(i, j, k1);
            }
    }
}

/// @brief 遍历一个BC设置引用的真实mesh set物理点。
template <typename T, typename Fn>
void forBoundarySettingPoints(Field& field,
                              const BCSetting<T>& bc,
                              Fn&& fn) {
    const auto& allSets = field.getAllSets();
    auto setIt = allSets.find(bc.name);
    if (setIt == allSets.end() || setIt->second.empty()) return;

    const int axis = (bc.type == EMPTY)
        ? Geometry::boundaryAxisForSet(field, bc.name)
        : Geometry::activeBoundaryAxisForSet(field, bc.name);

    const std::vector<int>& indices = setIt->second;
    for (int idx : indices) {
        int i = 0, j = 0, k = 0;
        field.getIJK(idx, i, j, k);
        if (!Geometry::isPhysicalPoint(field, i, j, k)) continue;
        fn(i, j, k, axis);
    }
}

/// @brief 应用固定值边界条件。
template <typename T>
void setFixedValue(Field& field, int i, int j, int k,
                   int axis, const T& bcValue, int vIdx,
                   bool useILW = false, int ilwOrder = 0) {
    const auto reconstruction =
        Reconstruction::fromILWSetting(useILW, ilwOrder);
    if constexpr (std::is_same_v<T, double>) {
        Reconstruction::applyScalar(
            field, i, j, k, axis, vIdx, bcValue,
            Law::Kind::FixedValue, reconstruction);
    } else if constexpr (std::is_same_v<T, Vector3>) {
        Reconstruction::applyVector(
            field, i, j, k, axis, vIdx, bcValue,
            Law::Kind::FixedValue, reconstruction);
    } else {
        static_assert(AlwaysFalse<T>::value, "Unsupported boundary value type");
    }
}

/// @brief 应用零梯度边界条件。
template <typename T>
void setZeroGradient(Field& field, int i, int j, int k,
                     int axis, int vIdx, bool useILW = false,
                     int ilwOrder = 0) {
    const auto reconstruction =
        Reconstruction::fromILWSetting(useILW, ilwOrder);
    if constexpr (std::is_same_v<T, double>) {
        Reconstruction::applyScalar(
            field, i, j, k, axis, vIdx, 0.0,
            Law::Kind::ZeroGradient, reconstruction);
    } else if constexpr (std::is_same_v<T, Vector3>) {
        Reconstruction::applyVector(
            field, i, j, k, axis, vIdx, Vector3(),
            Law::Kind::ZeroGradient, reconstruction);
    } else {
        static_assert(AlwaysFalse<T>::value, "Unsupported boundary value type");
    }
}

/// @brief 应用对称/滑移边界条件。
template <typename T>
void setSymmetry(Field& field, int i, int j, int k,
                 int axis, int vIdx, bool useILW = false,
                 int ilwOrder = 0) {
    const auto reconstruction =
        Reconstruction::fromILWSetting(useILW, ilwOrder);
    if constexpr (std::is_same_v<T, double>) {
        Reconstruction::applyScalar(
            field, i, j, k, axis, vIdx, 0.0,
            Law::Kind::Symmetry, reconstruction);
    } else if constexpr (std::is_same_v<T, Vector3>) {
        Reconstruction::applyVector(
            field, i, j, k, axis, vIdx, Vector3(),
            Law::Kind::Symmetry, reconstruction);
    } else {
        static_assert(AlwaysFalse<T>::value, "Unsupported boundary value type");
    }
}

/// @brief 应用empty边界条件。
template <typename T>
void setEmpty(Field& field, int i, int j, int k,
              int vIdx, int axis, bool useILW = false) {
    const auto reconstruction = useILW
        ? Reconstruction::Selection{Reconstruction::Kind::ILW, 0}
        : Reconstruction::Selection{Reconstruction::Kind::Linear, 0};
    if constexpr (std::is_same_v<T, double>) {
        Reconstruction::applyScalar(
            field, i, j, k, axis, vIdx, 0.0,
            Law::Kind::Empty, reconstruction);
    } else if constexpr (std::is_same_v<T, Vector3>) {
        Reconstruction::applyVector(
            field, i, j, k, axis, vIdx, Vector3(),
            Law::Kind::Empty, reconstruction);
    } else {
        static_assert(AlwaysFalse<T>::value, "Unsupported boundary value type");
    }
}

/// @brief 按配置集合更新一个变量族的边界条件。
template <typename T>
void update(Field& field,
            const std::vector<BCSetting<T>>& bcSettings,
            int vIdx,
            bool useILW = false,
            int ilwOrder = 0) {
    for (const auto& bc : bcSettings) {
        forBoundarySettingPoints(field, bc, [&](int i, int j, int k, int axis) {
            switch (bc.type) {
                case FIXED_VALUE:
                    setFixedValue<T>(field, i, j, k, axis, bc.value,
                                     vIdx, useILW, ilwOrder);
                    break;
                case ZERO_GRADIENT:
                    setZeroGradient<T>(field, i, j, k, axis, vIdx,
                                       useILW, ilwOrder);
                    break;
                case SYMMETRY:
                    setSymmetry<T>(field, i, j, k, axis, vIdx,
                                   useILW, ilwOrder);
                    break;
                case EMPTY:
                    setEmpty<T>(field, i, j, k, vIdx, axis, useILW);
                    break;
                default:
                    break;
            }
        });
    }
}

} // namespace Boundary
} // namespace SF
