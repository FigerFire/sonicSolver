#pragma once

/// @file SF_stateRealizer.h
/// @brief Resolved unknowns 到 StateBundle 已有 storage 的启动阶段绑定。

#include "SF_resolvedSimulationSystem.h"
#include "core/state/SF_stateBundle.h"

#include <vector>

namespace SF::System {

/// @brief 一个已解析未知量的稳定 runtime handle。
struct RealizedUnknown {
    const UnknownDescriptor* descriptor = nullptr;
    std::vector<State::DistributedFieldView*> fields;
};

/// @brief 已完成验证的非拥有 state contract。
class StateRealization {
public:
    const RealizedUnknown& at(std::string_view id) const;
    std::size_t size() const { return unknowns_.size(); }

private:
    friend StateRealization realizeState(
        const ResolvedSimulationSystem&, State::StateBundle&);
    std::vector<RealizedUnknown> unknowns_;
};

/// @brief 一次性解析并验证 unknown storage；不分配或复制数值数组。
StateRealization realizeState(
    const ResolvedSimulationSystem& system, State::StateBundle& state);

} // namespace SF::System
