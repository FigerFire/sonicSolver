#pragma once
#include <string>
#include <vector>
namespace SF::System {
// A version is an algorithm input, not an instruction to refresh a value.
enum class StateVersion { Current, Stage, OldTime, Frozen, Lagged };
enum class StateEvaluation { Direct, Lazy, Materialized };
struct StateUse {
    std::string symbol;
    StateVersion version=StateVersion::Current;
    bool halo=false;
    std::string snapshot;
    bool perIteration=false;
};
struct StateEffect {
    std::string symbol;
    bool publish=false;
    bool halo=false;
    bool invalidatesLazy=true;
};
/// @brief 能力必须服务于同一方程、流体端口和边界；空绑定仅表示全局能力。
struct CapabilityBinding {
    std::string name;
    std::string equation;
    std::string fluidPort;
    std::string boundary;
};
struct CapabilityRequirement {
    std::string when;
    std::string required;
    std::string reason;
    std::string equation;
    std::string fluidPort;
    std::string boundary;
};
// Edges validate a recipe's chosen order. They never select a coupling algorithm.
struct PlacementRequirement {
    std::string scope;
    std::string occurrence;
    std::vector<std::string> after;
    std::vector<std::string> before;
};
inline const char* toString(StateVersion version) {
    switch(version) {
    case StateVersion::Current:return "current";
    case StateVersion::Stage:return "stage";
    case StateVersion::OldTime:return "old-time";
    case StateVersion::Frozen:return "frozen";
    case StateVersion::Lagged:return "lagged";
    }
    return "unknown";
}
}
