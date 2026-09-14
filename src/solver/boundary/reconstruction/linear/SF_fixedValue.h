/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_fixedValue.h
/// @brief 代数型固定值边界条件。

#include "SF_field.h"

namespace SF {
namespace Boundary {
namespace Algebraic {
namespace FixedValue {

/// @brief 应用标量固定值边界。
void applyScalar(Field& field, int i, int j, int k,
                 int axis, double bcValue, int vIdx);

/// @brief 应用Vector3固定值边界。
void applyVector3(Field& field, int i, int j, int k,
                  int axis, const SF::Vector3& bcValue, int vIdx);

} // namespace FixedValue
} // namespace Algebraic
} // namespace Boundary
} // namespace SF
