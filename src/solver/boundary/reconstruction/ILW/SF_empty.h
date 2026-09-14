/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_empty.h
/// @brief ILW empty边界条件。

#include "SF_field.h"

namespace SF {
namespace Boundary {
namespace ILW {
namespace Empty {

/// @brief 应用ILW标量empty边界，指定轴法向导数为0。
void applyScalar(Field& field, int i, int j, int k, int vIdx, int axis);

/// @brief 应用ILW Vector3 empty边界，指定轴法向导数为0。
void applyVector3(Field& field, int i, int j, int k, int vIdx, int axis);

} // namespace Empty
} // namespace ILW
} // namespace Boundary
} // namespace SF
