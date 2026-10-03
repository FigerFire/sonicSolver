#pragma once

/// @file SF_termKernel.h
/// @brief Neutral compiled source operation on conservative state and residual.

#include <functional>
#include "core/state/SF_valueTypes.h"

namespace SF {
class Field;
class Residual;

namespace System {
using ConservativeSourceKernel = std::function<void(Field&,Residual&)>;
using PrimitiveMomentumSource = std::function<Vector3(
    const Field&,int,int,int,const Vector3&)>;
}
} // namespace SF
