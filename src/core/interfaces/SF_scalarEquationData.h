#pragma once
#include <array>
#include <functional>
#include <map>
#include <string>
#include <vector>
#include <memory>

namespace SF::FDM {
using ScalarKnownValue = std::function<double(double,double,double,double)>;
enum class ScalarBoundaryKind { Unspecified, FixedValue, ZeroGradient };
struct ScalarBoundaryCondition {
    ScalarBoundaryKind kind=ScalarBoundaryKind::Unspecified;
    ScalarKnownValue value;
};
/// Non-owning numerical data injection. Storage/clock remain in StateBundle.
/// Face order: x-/x+, y-/y+, z-/z+. Known functions receive x,y,z,stageTime.
struct ScalarEquationData {
    std::string target;
    std::array<ScalarBoundaryCondition,6> boundary;
    std::map<std::string,ScalarKnownValue> sources;
};
/// @brief Caller keeps each immutable data object alive throughout execution.
struct ScalarInstanceData {
    std::string occurrence;
    std::weak_ptr<const ScalarEquationData> data;
};
}
