/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_zeroGradient.h
/// @brief 高阶ILW零梯度边界条件。
///
/// 基于多维张量积 WENO 型高阶逆 Lax-Wendroff 方法：
///   - 多维模板外推：3/5/7/9 阶模板 + 非线性光滑权重
///   - 零梯度约束 ∂U/∂n = 0 直接消去一阶法向导数
///   - 二阶及以上导数由 WENO 外推给出
///   - Taylor 展开构造各层虚胞值
///
/// `[numerics].ILW=3/5/7/9` 分别使用 2/4/6/8 阶 Taylor 展开。
/// 样本不足时终止并报告位置，不做静默降阶兜底。

#include "SF_field.h"

namespace SF {
namespace Boundary {
namespace ILW {
namespace ZeroGradient {

/// @brief 应用ILW标量零梯度边界。
void applyScalar(Field& field, int i, int j, int k, int vIdx,
                 int accuracyOrder);

/// @brief 应用ILW Vector3零梯度边界。
void applyVector3(Field& field, int i, int j, int k, int vIdx,
                  int accuracyOrder);

} // namespace ZeroGradient
} // namespace ILW
} // namespace Boundary
} // namespace SF
