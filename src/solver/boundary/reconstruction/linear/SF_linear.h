/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_linear.h
/// @brief 线性/代数型物理边界重构调度入口。

#include "SF_field.h"

namespace SF {
namespace Boundary {
namespace Algebraic {

/// @brief 固定值标量边界调度。
void setFixedValueScalar(Field& field, int i, int j, int k,
                         int axis, double bcValue, int vIdx);

/// @brief 固定值Vector3边界调度。
void setFixedValueVector3(Field& field, int i, int j, int k,
                          int axis, const SF::Vector3& bcValue, int vIdx);

/// @brief 零梯度标量边界调度。
void setZeroGradientScalar(Field& field, int i, int j, int k, int axis, int vIdx);

/// @brief 零梯度Vector3边界调度。
void setZeroGradientVector3(Field& field, int i, int j, int k, int axis, int vIdx);

/// @brief 对称标量边界调度。
void setSymmetryScalar(Field& field, int i, int j, int k, int axis, int vIdx);

/// @brief 对称Vector3边界调度。
void setSymmetryVector3(Field& field, int i, int j, int k, int axis, int vIdx);

/// @brief empty标量边界调度。
void setEmptyScalar(Field& field, int i, int j, int k, int vIdx, int axis);

/// @brief empty Vector3边界调度。
void setEmptyVector3(Field& field, int i, int j, int k, int vIdx, int axis);

} // namespace Algebraic
} // namespace Boundary
} // namespace SF
