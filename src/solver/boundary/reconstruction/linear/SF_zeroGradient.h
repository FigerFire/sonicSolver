/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_zeroGradient.h
/// @brief 代数型零梯度边界条件。

#include "SF_field.h"

namespace SF {
namespace Boundary {
namespace Algebraic {
namespace ZeroGradient {

/// @brief 应用标量零梯度边界。
void applyScalar(Field& field, int i, int j, int k, int axis, int vIdx);

/// @brief 应用Vector3零梯度边界。
void applyVector3(Field& field, int i, int j, int k, int axis, int vIdx);

} // namespace ZeroGradient
} // namespace Algebraic
} // namespace Boundary
} // namespace SF
