#pragma once

/// @file SF_assemblyPlan.h
/// @brief Equation::Definition 到 specialized assembler 的非拥有 lowering。

#include "core/system/SF_expression.h"

#include <string>
#include <string_view>
#include <vector>

namespace SF::Equation {

/// @brief 单个 equation 的 executable operator view；不复制 physics metadata。
struct AssemblyPlan {
    const Definition* definition = nullptr;
    std::vector<const Term*> left;
    std::vector<const Term*> right;

    bool contains(TermKind kind) const;
};

/// @brief 启动阶段一次性建立的 equation-id -> AssemblyPlan 稳定表。
class AssemblyPlanRegistry {
public:
    AssemblyPlanRegistry() = default;
    AssemblyPlanRegistry(
        const System& definitions,
        const std::vector<std::string>& equationIds);

    const AssemblyPlan& at(std::string_view equationId) const;
    bool contains(std::string_view equationId) const;
    std::size_t size() const { return plans_.size(); }

private:
    struct Entry {
        std::string id;
        AssemblyPlan plan;
    };
    std::vector<Entry> plans_;
};

/// @brief 从 authoritative definitions 建立稳定的非拥有 assembly views。
std::vector<AssemblyPlan> makeAssemblyPlan(
    const System& definitions,
    const std::vector<std::string>& equationIds);

} // namespace SF::Equation
