#pragma once

/// @file SF_termKernel.h
/// @brief Neutral compiled source operation on conservative state and residual.

#include <functional>
#include <vector>
#include <string>
#include "core/state/SF_valueTypes.h"

namespace SF {
class Field;
class Residual;

namespace System {
struct StageConservativeSource {
    using Reader=std::function<double(int,int)>;
    std::vector<std::string> reads;
    std::function<void(Field&,Residual&,const std::vector<Reader>&)> evaluate;
};
using ConservativeSourceKernel = std::function<void(Field&,Residual&)>;
using PrimitiveMomentumSource = std::function<Vector3(
    const Field&,int,int,int,const Vector3&)>;
}
} // namespace SF
