#pragma once

/// @file SF_termKernel.h
/// @brief Neutral compiled source operation on conservative state and residual.

#include <functional>
#include <vector>
#include <string>
#include "core/state/SF_valueTypes.h"
#include "core/system/SF_dataContract.h"

namespace SF {
class Field;
class Residual;

namespace System {
/// @brief WHAT source(expr) 的冻结求值；绑定 STATE reader，不拥有更新或时间循环。
struct FieldSource {
    using Reader=std::function<double(int,int)>;
    std::vector<StateUse> inputs;
    std::vector<std::string> evaluation;
    std::function<Reader(std::vector<Reader>)> bind;
};
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
