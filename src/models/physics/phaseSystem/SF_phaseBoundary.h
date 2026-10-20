#pragma once

/// @file SF_phaseBoundary.h
/// @brief 独立相标量和速度边界条件实现。

#include "SF_phaseProperties.h"
#include "SF_phaseState.h"

namespace SF::Physics::PhaseSystems {

/// @brief 每相原始量边界与主状态提交的独立应用器。
class PhaseBoundaryApplicator {
public:
    /// @brief 将相初值写入 primitive cache，随后构造主状态。
    void initialize(const Field& geometry,
                    const Multiphase::PhaseProperties& properties,
                    double initialPressure,
                    PhaseState& state) const;

    /// @brief 应用相原始量物理边界和 ghost 延拓，并提交主状态。
    void apply(const Field& geometry,
               const Multiphase::PhaseProperties& properties,
               PhaseState& state) const;

    /// @brief 应用共享压力边界并闭合 ghost；不改写任一相的主状态。
    void applySharedPressure(
        const Field& geometry,
        ScalarField& pressure,
        const std::vector<BCSetting<double>>& settings) const;
};

} // namespace SF::Physics::PhaseSystems
