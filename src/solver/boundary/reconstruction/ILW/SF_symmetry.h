/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_symmetry.h
/// @brief ILW对称/滑移边界条件。

#include "SF_field.h"

namespace SF {
namespace Boundary {
namespace ILW {
namespace Symmetry {

/// @brief 应用ILW标量对称边界。
void applyScalar(Field& field, int i, int j, int k, int vIdx,
                 int accuracyOrder);

/// @brief 应用ILW Vector3对称边界。
void applyVector3(Field& field, int i, int j, int k, int vIdx,
                  int accuracyOrder);

} // namespace Symmetry
} // namespace ILW
} // namespace Boundary
} // namespace SF
