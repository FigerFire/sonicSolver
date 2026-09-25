/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_cellType.h
/// @brief IBM单元类型判定工具。
///
/// 本文件只集中维护IBM相关的CellFlag语义，不做通量、边界或ILW算法。

#include "SF_field.h"

namespace SF {
namespace IBM {

/// @brief 判断索引是否位于物理网格区域内。
/// @param field 结构网格场。
/// @param i i方向全局存储索引。
/// @param j j方向全局存储索引。
/// @param k k方向全局存储索引。
/// @return 位于不含外层ghost的物理区域时返回true。
inline bool isPhysicalCell(const Field& field, int i, int j, int k) {
    const int ng = field.NG();
    return i >= ng && i < ng + field.NX()
        && j >= ng && j < ng + field.NY()
        && k >= ng && k < ng + field.NZ();
}

/// @brief 判断单元数据是否可供IBM/ILW局部闭合读取。
/// @param field 结构网格场。
/// @param i i方向全局存储索引。
/// @param j j方向全局存储索引。
/// @param k k方向全局存储索引。
/// @return 物理点或已建立通信映射的MPI halo点返回true。
inline bool isIBMReadableCell(const Field& field, int i, int j, int k) {
    if (i < 0 || i >= field.MX()
        || j < 0 || j >= field.MY()
        || k < 0 || k >= field.MZ()) {
        return false;
    }
    return isPhysicalCell(field, i, j, k)
        || field.isCommunicationHalo(i, j, k);
}

/// @brief 判断单元是否是IBM内部的非流体单元。
/// @param field 结构网格场。
/// @param i i方向全局存储索引。
/// @param j j方向全局存储索引。
/// @param k k方向全局存储索引。
/// @return 物理区域内的IBM_GHOST_CELL或SOLID_CELL返回true。
inline bool isIbmNonFluidCell(const Field& field, int i, int j, int k) {
    return isPhysicalCell(field, i, j, k)
        && field.CellFlag(i, j, k) != FLUID_CELL;
}

/// @brief 判断单元是否是物理区域内的真实流体单元。
/// @param field 结构网格场。
/// @param i i方向全局存储索引。
/// @param j j方向全局存储索引。
/// @param k k方向全局存储索引。
/// @return 物理区域内且CellFlag为FLUID_CELL时返回true。
inline bool isFluidCell(const Field& field, int i, int j, int k) {
    return isPhysicalCell(field, i, j, k)
        && field.CellFlag(i, j, k) == FLUID_CELL;
}

} // namespace IBM
} // namespace SF
