/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_empty.h
/// @brief 代数型empty边界条件。

#include "SF_field.h"

namespace SF {
namespace Boundary {
namespace Algebraic {
namespace Empty {

/// @brief 应用标量empty边界，沿指定轴采用常数外推。
void applyScalar(Field& field, int i, int j, int k, int vIdx, int axis);

/// @brief 应用Vector3 empty边界，沿指定轴采用常数外推。
void applyVector3(Field& field, int i, int j, int k, int vIdx, int axis);

} // namespace Empty
} // namespace Algebraic
} // namespace Boundary
} // namespace SF
