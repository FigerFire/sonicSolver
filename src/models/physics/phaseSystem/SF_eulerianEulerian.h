#pragma once

/// @file SF_eulerianEulerian.h
/// @brief 每相独立 U/rho/T、共享压力的 Eulerian–Eulerian PhaseSystem。

#include "SF_phaseSystem.h"
#include "SF_phaseBoundary.h"

namespace SF::Physics::PhaseSystems {

class EulerianEulerianPhaseSystem final : public PhaseSystem {
public:
    explicit EulerianEulerianPhaseSystem(
        Multiphase::MultiPhaseConfig config);

    void initialize(const Field& geometry) override;
    void recoverPrimitiveState() override;
    void reconstructReferencePhaseMass() override;
    void computeInterphase(
        double dt,
        const std::vector<PhaseVectorField>& previousVelocity) override;
    void validateState(const std::string& stage) const override;

    const Field& geometry() const override;
    ScalarField& sharedPressure() override { return pressure_; }
    const ScalarField& sharedPressure() const override { return pressure_; }
    std::vector<PhaseState>& phases() override { return phases_; }
    const std::vector<PhaseState>& phases() const override { return phases_; }
    PhaseEquationSources& sources() override { return sources_; }
    const PhaseEquationSources& sources() const override { return sources_; }
    size_t referencePhaseIndex() const override { return referencePhase_; }
    const Multiphase::PhaseProperties& phaseProperties(
        size_t phase) const override;
    const Multiphase::MultiPhaseConfig& config() const override { return config_; }

private:
    Multiphase::MultiPhaseConfig config_;
    const Field* geometry_ = nullptr;
    std::vector<PhaseState> phases_;
    std::vector<const Multiphase::PhaseProperties*> properties_;
    ScalarField pressure_;
    PhaseEquationSources sources_;
    size_t referencePhase_ = 0;
    std::vector<std::array<size_t, 2>> phasePairs_;
    PhaseBoundaryApplicator boundaryApplicator_;

    size_t phaseIndex(const std::string& name) const;
    void recoverPhase(size_t phase);
};

} // namespace SF::Physics::PhaseSystems
