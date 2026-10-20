#pragma once

/// @file SF_stateLayout.h
/// @brief 方程组主变量与算法更新策略的值类型描述。

#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace SF::State {

/// @brief 变量所在的离散位置。
enum class FieldLocation { Cell, Face, Node };

/// @brief 变量的守恒语义。
enum class ConservationKind {
    Conservative,
    BoundedConservative,
    AdvectedGeometry,
    Algorithmic
};

/// @brief workflow 对变量采用的更新方式。
enum class UpdatePolicy {
    Explicit,
    BoundedExplicit,
    ExternalExplicit,
    HamiltonJacobi,
    LocalImplicit,
    PressureCorrection,
    DerivedOnly
};

/// @brief transported state 需要 halo 同步的离散阶段。
enum class HaloSyncStage : unsigned {
    None = 0,
    BeforeRHS = 1u << 0,
    AfterUpdate = 1u << 1,
    BeforeOutput = 1u << 2
};

inline HaloSyncStage operator|(
        HaloSyncStage first, HaloSyncStage second) {
    return static_cast<HaloSyncStage>(
        static_cast<unsigned>(first)
        | static_cast<unsigned>(second));
}

inline bool contains(
        HaloSyncStage stages, HaloSyncStage stage) {
    return (static_cast<unsigned>(stages)
            & static_cast<unsigned>(stage)) != 0;
}

/// @brief 一个变量组的 halo 深度和同步频率契约。
struct HaloPolicy {
    int depth = 0;
    HaloSyncStage stages = HaloSyncStage::None;

    void validate(const std::string& variable) const {
        if (depth < 0) {
            throw std::runtime_error(
                "StateLayout variable '" + variable
                + "' has a negative halo depth.");
        }
        if (depth == 0 && stages != HaloSyncStage::None) {
            throw std::runtime_error(
                "StateLayout variable '" + variable
                + "' requests synchronization with zero halo depth.");
        }
    }
};

/// @brief 可选的变量范围；disabled 时上下界不参与更新。
struct Bounds {
    bool enabled = false;
    double lower = 0.0;
    double upper = 0.0;

    /// @brief 检查边界值本身是否合法，不修改变量。
    void validate(const std::string& variable) const {
        if (enabled && (!std::isfinite(lower) || !std::isfinite(upper)
                        || lower > upper)) {
            throw std::runtime_error(
                "StateLayout variable '" + variable
                + "' has invalid ordered finite bounds.");
        }
    }
};

/// @brief 一个主变量或 transported variable 的静态契约。
struct VariableDescriptor {
    std::string name;
    FieldLocation location = FieldLocation::Cell;
    ConservationKind conservation = ConservationKind::Conservative;
    UpdatePolicy updatePolicy = UpdatePolicy::Explicit;
    Bounds bounds;
    std::string unit;
    HaloPolicy halo;
};

/// @brief FluidStateModel 提供给 workflow 的变量布局。
class StateLayout {
public:
    StateLayout() = default;
    explicit StateLayout(std::vector<VariableDescriptor> variables)
        : variables_(std::move(variables)) {
        validate();
    }

    /// @brief 返回全部变量描述。
    const std::vector<VariableDescriptor>& variables() const {
        return variables_;
    }

    /// @brief 返回布局中的变量数。
    int size() const { return static_cast<int>(variables_.size()); }

    /// @brief 按名称查找变量；不存在返回 -1。
    int find(const std::string& name) const {
        for (int i = 0; i < size(); ++i) {
            if (variables_[(size_t)i].name == name) return i;
        }
        return -1;
    }

    /// @brief 检查空名称、重复名称和非法 bounds。
    void validate() const {
        for (size_t i = 0; i < variables_.size(); ++i) {
            const auto& variable = variables_[i];
            if (variable.name.empty()) {
                throw std::runtime_error("StateLayout contains an unnamed variable.");
            }
            variable.bounds.validate(variable.name);
            variable.halo.validate(variable.name);
            for (size_t j = 0; j < i; ++j) {
                if (variables_[j].name == variable.name) {
                    throw std::runtime_error(
                        "StateLayout contains duplicate variable '"
                        + variable.name + "'.");
                }
            }
        }
    }

private:
    std::vector<VariableDescriptor> variables_;
};

} // namespace SF::State
