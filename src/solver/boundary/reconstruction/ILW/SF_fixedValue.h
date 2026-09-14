/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_fixedValue.h
/// @brief 高阶 ILW 固定值边界条件。
///
/// 基于多维张量积 WENO 型边界外推：
///   - U^(0) = 预设边界值
///   - 多维模板外推计算高阶法向导数
///   - Taylor 展开构造各层虚胞值
///
/// `[numerics].ILW=3/5/7/9` 分别使用 2/4/6/8 阶 Taylor 展开。
/// 样本不足时终止并报告位置，不做静默降阶兜底。

#include "SF_field.h"

namespace SF {
namespace Boundary {
namespace ILW {
namespace FixedValue {

/// @brief 应用ILW标量固定值边界。
void applyScalar(Field& field, int i, int j, int k,
                 double bcValue, int vIdx, int accuracyOrder);

/// @brief 应用ILW Vector3固定值边界。
void applyVector3(Field& field, int i, int j, int k,
                  const SF::Vector3& bcValue, int vIdx,
                  int accuracyOrder);

} // namespace FixedValue
} // namespace ILW
} // namespace Boundary
} // namespace SF
