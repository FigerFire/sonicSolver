#pragma once

/// @file SF_variableRegistry.h
/// @brief 与主守恒状态同步推进的通用标量变量注册表。

#include "SF_scalarField.h"
#include "SF_distributedField.h"
#include "SF_stateLayout.h"

#include <stdexcept>
#include <utility>
#include <vector>

namespace SF::State {

/// @brief 一个由时间积分器消费的 transported scalar 视图。
struct IntegratedVariable {
    VariableDescriptor descriptor;
    ScalarField* value = nullptr;
    std::vector<double>* rhs = nullptr;
};

/// @brief 管理任意数量 transported scalar 的非拥有注册表。
class VariableRegistry {
public:
    /// @brief 注册一个值场及其 RHS；注册表不拥有二者。
    void add(VariableDescriptor descriptor,
             ScalarField& value,
             std::vector<double>& rhs) {
        descriptor.bounds.validate(descriptor.name);
        if (descriptor.name.empty()) {
            throw std::runtime_error("VariableRegistry cannot register an unnamed variable.");
        }
        for (const auto& item : variables_) {
            if (item.descriptor.name == descriptor.name) {
                throw std::runtime_error(
                    "VariableRegistry duplicate variable '" + descriptor.name + "'.");
            }
        }
        variables_.push_back({std::move(descriptor), &value, &rhs});
    }

    /// @brief 注册由专用算法推进或纯派生的标量；不要求通用时间积分 RHS。
    void add(VariableDescriptor descriptor,
             ScalarField& value) {
        if (descriptor.updatePolicy == UpdatePolicy::Explicit
            || descriptor.updatePolicy == UpdatePolicy::BoundedExplicit
            || descriptor.updatePolicy == UpdatePolicy::HamiltonJacobi) {
            throw std::runtime_error(
                "VariableRegistry tableau-integrated variable '" + descriptor.name
                + "' requires an RHS view.");
        }
        descriptor.bounds.validate(descriptor.name);
        if (descriptor.name.empty()) {
            throw std::runtime_error(
                "VariableRegistry cannot register an unnamed variable.");
        }
        for (const auto& item : variables_) {
            if (item.descriptor.name == descriptor.name) {
                throw std::runtime_error(
                    "VariableRegistry duplicate variable '"
                    + descriptor.name + "'.");
            }
        }
        variables_.push_back({std::move(descriptor), &value, nullptr});
    }

    /// @brief 清空非拥有视图，不修改任何场数据。
    void clear() { variables_.clear(); }

    bool empty() const { return variables_.empty(); }
    size_t size() const { return variables_.size(); }

    const std::vector<IntegratedVariable>& variables() const {
        return variables_;
    }

    /// @brief 将所有已登记辅助量发布到统一分布式状态注册表。
    void registerDistributed(
            DistributedFieldRegistry& distributed,
            int blockId, Field& geometry) const {
        for (const auto& variable : variables_) {
            distributed.add(scalarView(
                variable.descriptor.name, blockId, geometry,
                *variable.value));
        }
    }

    /// @brief 检查所有 value/RHS 尺寸与主 Field 一致。
    void validateFor(const Field& field) const {
        for (const auto& item : variables_) {
            if (!item.value) {
                throw std::runtime_error(
                    "VariableRegistry contains a null value view.");
            }
            if (!item.value->isCompatibleWith(field)) {
                throw std::runtime_error(
                    "Transported variable '" + item.descriptor.name
                    + "' dimensions do not match Field.");
            }
            const UpdatePolicy policy = item.descriptor.updatePolicy;
            const bool tableauIntegrated =
                policy == UpdatePolicy::Explicit
                || policy == UpdatePolicy::BoundedExplicit
                || policy == UpdatePolicy::HamiltonJacobi;
            if (tableauIntegrated && !item.rhs) {
                throw std::runtime_error(
                    "Transported variable '" + item.descriptor.name
                    + "' requires an RHS for unified explicit-tableau integration.");
            }
            if (item.rhs
                && item.rhs->size() != item.value->values().size()) {
                throw std::runtime_error(
                    "Transported variable '" + item.descriptor.name
                    + "' RHS size does not match its value field.");
            }
            if (policy != UpdatePolicy::Explicit
                && policy != UpdatePolicy::BoundedExplicit
                && policy != UpdatePolicy::ExternalExplicit
                && policy != UpdatePolicy::HamiltonJacobi
                && policy != UpdatePolicy::LocalImplicit
                && policy != UpdatePolicy::PressureCorrection
                && policy != UpdatePolicy::DerivedOnly) {
                throw std::runtime_error(
                    "VariableRegistry received an unsupported update policy for '"
                    + item.descriptor.name + "'.");
            }
        }
    }

private:
    std::vector<IntegratedVariable> variables_;
};

} // namespace SF::State
