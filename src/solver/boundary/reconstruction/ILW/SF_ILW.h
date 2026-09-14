/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_ILW.h
/// @brief ILW物理边界调度入口。

#include "SF_field.h"
#include "core/mesh/SF_dimension.h"
#include "SF_cellType.h"

#include <algorithm>

namespace SF {
namespace Boundary {
namespace ILW {

// ═══════════════════════════════════════════════════════════════
//  ILW 物理边界共享工具
// ═══════════════════════════════════════════════════════════════

/// @brief 判断索引是否位于物理域外的普通虚胞层。
inline bool isGhostCell(const Field& field, int i, int j, int k) {
    const int ng = field.NG();
    return (i < ng || i >= field.NX() + ng ||
            j < ng || j >= field.NY() + ng ||
            k < ng || k >= field.NZ() + ng);
}

/// @brief 将索引夹到物理域内最近的实胞。
inline void nearestInteriorCell(const Field& field,
                                int i, int j, int k,
                                int& ri, int& rj, int& rk) {
    const int ng = field.NG();
    ri = std::max(ng, std::min(i, field.NX() + ng - 1));
    rj = std::max(ng, std::min(j, field.NY() + ng - 1));
    rk = std::max(ng, std::min(k, field.NZ() + ng - 1));
}

/// @brief 返回普通虚胞点的主法向轴，0/1/2对应x/y/z，-1表示非虚胞。
inline int boundaryAxis(const Field& field, int i, int j, int k) {
    const int ng = field.NG();
    if (i < ng || i >= field.NX() + ng) return 0;
    if (j < ng || j >= field.NY() + ng) return 1;
    if (k < ng || k >= field.NZ() + ng) return 2;
    return -1;
}

/// @brief 反射速度/动量的法向分量。
inline Vector3 reflectNormalComponent(const Field& field,
                                      int i, int j, int k,
                                      const Vector3& value) {
    Vector3 reflected = value;
    const int axis = boundaryAxis(field, i, j, k);
    if (axis == 0) reflected.x = -reflected.x;
    if (axis == 1) reflected.y = -reflected.y;
    if (axis == 2) reflected.z = -reflected.z;
    return reflected;
}

/// @brief 对实边界点外侧的普通虚胞层执行广播。
///
/// 只填充活跃方向的外层虚胞；被 EMPTY 关闭的方向不在此处填充。
template <typename Fn>
void forBoundaryGhosts(Field& field, int i, int j, int k, Fn&& fn) {
    const int ng = field.NG();
    const int nx = field.NX();
    const int ny = field.NY();
    const int nz = field.NZ();

    if (Math::isDirectionActiveIndex(0)) {
        if (i == ng) {
            for (int b = 0; b < ng; ++b) fn(b, j, k, i, j, k);
        }
        if (i == nx + ng - 1) {
            for (int b = 1; b <= ng; ++b) fn(i + b, j, k, i, j, k);
        }
    }
    if (Math::isDirectionActiveIndex(1)) {
        if (j == ng) {
            for (int b = 0; b < ng; ++b) fn(i, b, k, i, j, k);
        }
        if (j == ny + ng - 1) {
            for (int b = 1; b <= ng; ++b) fn(i, j + b, k, i, j, k);
        }
    }
    if (Math::isDirectionActiveIndex(2)) {
        if (k == ng) {
            for (int b = 0; b < ng; ++b) fn(i, j, b, i, j, k);
        }
        if (k == nz + ng - 1) {
            for (int b = 1; b <= ng; ++b) fn(i, j, k + b, i, j, k);
        }
    }
}

/// @brief 只沿指定普通物理边界法向填充外侧虚胞层。
///
/// 用于按当前 patch 法向应用 ILW 边界条件，避免角点处跨 patch 写入 ghost。
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

// ═══════════════════════════════════════════════════════════════
//  物理边界条件调度
// ═══════════════════════════════════════════════════════════════

/// @brief 固定值标量边界调度。
void setFixedValueScalar(Field& field, int i, int j, int k,
                         double bcValue, int vIdx, int accuracyOrder);

/// @brief 固定值Vector3边界调度。
void setFixedValueVector3(Field& field, int i, int j, int k,
                          const SF::Vector3& bcValue, int vIdx,
                          int accuracyOrder);

/// @brief 零梯度标量边界调度。
void setZeroGradientScalar(Field& field, int i, int j, int k, int vIdx,
                           int accuracyOrder);

/// @brief 零梯度Vector3边界调度。
void setZeroGradientVector3(Field& field, int i, int j, int k, int vIdx,
                            int accuracyOrder);

/// @brief 对称标量边界调度。
void setSymmetryScalar(Field& field, int i, int j, int k, int vIdx,
                       int accuracyOrder);

/// @brief 对称Vector3边界调度。
void setSymmetryVector3(Field& field, int i, int j, int k, int vIdx,
                        int accuracyOrder);

/// @brief empty标量边界调度。
void setEmptyScalar(Field& field, int i, int j, int k, int vIdx, int axis);

/// @brief empty Vector3边界调度。
void setEmptyVector3(Field& field, int i, int j, int k, int vIdx, int axis);

} // namespace ILW

/// @brief 判断索引是否位于物理网格区域内（不含ghost层）。
inline bool isPhysicalCell(const Field& field, int i, int j, int k) {
    return IBM::isPhysicalCell(field, i, j, k);
}

/// @brief 判断单元是否是IBM内部的非流体单元。
inline bool isIbmNonFluidCell(const Field& field, int i, int j, int k) {
    return IBM::isIbmNonFluidCell(field, i, j, k);
}

/// @brief 从Field拷贝一个守恒状态。
inline void loadConservative(const Field& field, int i, int j, int k, double q[5]) {
    q[0] = field(i, j, k, RHO);
    q[1] = field(i, j, k, RU);
    q[2] = field(i, j, k, RV);
    q[3] = field(i, j, k, RW);
    q[4] = field(i, j, k, E);
}

} // namespace Boundary
} // namespace SF
