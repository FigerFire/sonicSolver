/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_symmetry.h
/// @brief 代数型对称边界条件。

#include "SF_field.h"

namespace SF {
namespace Boundary {
namespace Algebraic {
namespace Symmetry {

/// @brief 应用标量对称边界。
void applyScalar(Field& field, int i, int j, int k, int axis, int vIdx);

/// @brief 应用Vector3对称边界。
void applyVector3(Field& field, int i, int j, int k, int axis, int vIdx);

} // namespace Symmetry
} // namespace Algebraic
} // namespace Boundary
} // namespace SF
