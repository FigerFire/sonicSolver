#pragma once

/// @file SF_capability.h
/// @brief 可组合求解模块使用的能力标识和值集合。

#include <initializer_list>
#include <set>
#include <string>

namespace SF::FDM {

/// @brief 一个模块能够提供或要求的数值/状态契约。
enum class Capability {
    CanonicalConservativeState,
    PrimitiveState,
    FiveVariableEulerState,
    CharacteristicEigenSystem,
    PhaseCanonicalState,
    ConservativeGhostState,
    PhaseGhostState,
    FittedBoundaryGeometry,
    ImmersedBoundaryGeometry,
    PredictedConservativeState,
    ImmersedConstraintProjection,
    LagrangeMultiplierField,
    AlgebraicBoundaryReconstruction,
    ILWBoundaryReconstruction,
    CanonicalFaceFlux,
    PressureCorrection,
    PressureJumpConsumer,
    TwoMaterialRiemann,
    SingleFieldExecution,
    MultiFieldExecution,
    FluidStateModelAwareTurbulence
};

/// @brief 返回稳定的能力名称，用于配置诊断和架构日志。
inline const char* toString(Capability capability) {
    switch (capability) {
        case Capability::CanonicalConservativeState:
            return "CanonicalConservativeState";
        case Capability::PrimitiveState: return "PrimitiveState";
        case Capability::FiveVariableEulerState:
            return "FiveVariableEulerState";
        case Capability::CharacteristicEigenSystem:
            return "CharacteristicEigenSystem";
        case Capability::PhaseCanonicalState: return "PhaseCanonicalState";
        case Capability::ConservativeGhostState:
            return "ConservativeGhostState";
        case Capability::PhaseGhostState: return "PhaseGhostState";
        case Capability::FittedBoundaryGeometry:
            return "FittedBoundaryGeometry";
        case Capability::ImmersedBoundaryGeometry:
            return "ImmersedBoundaryGeometry";
        case Capability::PredictedConservativeState:
            return "PredictedConservativeState";
        case Capability::ImmersedConstraintProjection:
            return "ImmersedConstraintProjection";
        case Capability::LagrangeMultiplierField:
            return "LagrangeMultiplierField";
        case Capability::AlgebraicBoundaryReconstruction:
            return "AlgebraicBoundaryReconstruction";
        case Capability::ILWBoundaryReconstruction:
            return "ILWBoundaryReconstruction";
        case Capability::CanonicalFaceFlux: return "CanonicalFaceFlux";
        case Capability::PressureCorrection: return "PressureCorrection";
        case Capability::PressureJumpConsumer:
            return "PressureJumpConsumer";
        case Capability::TwoMaterialRiemann: return "TwoMaterialRiemann";
        case Capability::SingleFieldExecution:
            return "SingleFieldExecution";
        case Capability::MultiFieldExecution:
            return "MultiFieldExecution";
        case Capability::FluidStateModelAwareTurbulence:
            return "FluidStateModelAwareTurbulence";
    }
    return "UnknownCapability";
}

/// @brief 去重的模块能力集合。
class CapabilitySet {
public:
    CapabilitySet() = default;
    CapabilitySet(std::initializer_list<Capability> values)
        : values_(values) {}

    void add(Capability capability) { values_.insert(capability); }
    void merge(const CapabilitySet& other) {
        values_.insert(other.values_.begin(), other.values_.end());
    }
    bool contains(Capability capability) const {
        return values_.find(capability) != values_.end();
    }
    bool empty() const { return values_.empty(); }
    const std::set<Capability>& values() const { return values_; }

    /// @brief 返回逗号分隔的稳定能力列表。
    std::string describe() const;

private:
    std::set<Capability> values_;
};

inline std::string CapabilitySet::describe() const {
    std::string result;
    for (Capability capability : values_) {
        if (!result.empty()) result += ", ";
        result += toString(capability);
    }
    return result;
}

} // namespace SF::FDM
