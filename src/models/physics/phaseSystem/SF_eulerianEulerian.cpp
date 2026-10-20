/// @file SF_eulerianEulerian.cpp
/// @brief PhaseSystem 状态、边界、源项与双欧拉物理实现。

#include "SF_eulerianEulerian.h"

#include "SF_interphase.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace SF::Physics::PhaseSystems {

EulerianEulerianPhaseSystem::EulerianEulerianPhaseSystem(
        Multiphase::MultiPhaseConfig config)
    : config_(std::move(config)) {
    Multiphase::validateMultiPhaseConfig(
        config_, "EulerianEulerianPhaseSystem");
}

size_t EulerianEulerianPhaseSystem::phaseIndex(
        const std::string& name) const {
    for (size_t phase = 0; phase < phases_.size(); ++phase) {
        if (Multiphase::normalizePhaseName(phases_[phase].name)
            == Multiphase::normalizePhaseName(name)) {
            return phase;
        }
    }
    throw std::runtime_error(
        "Eulerian phase name '" + name + "' is not declared.");
}

const Multiphase::PhaseProperties&
EulerianEulerianPhaseSystem::phaseProperties(size_t phase) const {
    if (phase >= properties_.size() || !properties_[phase]) {
        throw std::runtime_error(
            "Eulerian phase properties index is invalid.");
    }
    return *properties_[phase];
}

void EulerianEulerianPhaseSystem::initialize(const Field& geometry) {
    geometry_ = &geometry;
    phases_.clear();
    properties_.clear();
    pressure_.setupLike(geometry, "p", config_.initialPressure);

    for (const std::string& name :
         config_.eulerianEulerian.phaseNames) {
        const auto* properties = Multiphase::findPhase(config_, name);
        if (!properties) {
            throw std::runtime_error(
                "Eulerian phase properties are missing for '" + name + "'.");
        }
        phases_.emplace_back();
        phases_.back().setupLike(geometry, name);
        properties_.push_back(properties);
        boundaryApplicator_.initialize(
            geometry, *properties, config_.initialPressure, phases_.back());
    }

    referencePhase_ = phaseIndex(
        config_.eulerianEulerian.referencePhase);
    phasePairs_.clear();
    for (const auto& pair : config_.eulerianEulerian.phasePairs) {
        phasePairs_.push_back({
            phaseIndex(pair.continuousPhase),
            phaseIndex(pair.dispersedPhase)});
    }
    sources_.setupLike(geometry, phases_, phasePairs_);

    recoverPrimitiveState();
    for (size_t phase = 0; phase < phases_.size(); ++phase) {
        boundaryApplicator_.apply(
            geometry, phaseProperties(phase), phases_[phase]);
    }
    recoverPrimitiveState();
    validateState("initialize");
}

const Field& EulerianEulerianPhaseSystem::geometry() const {
    if (!geometry_) {
        throw std::runtime_error(
            "Eulerian PhaseSystem is not initialized.");
    }
    return *geometry_;
}

void EulerianEulerianPhaseSystem::recoverPhase(size_t phaseIndex) {
    PhaseState& phase = phases_.at(phaseIndex);
    const auto& properties = phaseProperties(phaseIndex);
    for (int cell = 0; cell < geometry().TotalSize(); ++cell) {
        const double mass =
            phase.primary.phaseMass.values()[(size_t)cell];
        const double pressure = pressure_.values()[(size_t)cell];
        const double phaseEnthalpy =
            phase.primary.phaseEnthalpy.values()[(size_t)cell];
        if (!std::isfinite(mass) || mass <= 0.0
            || !std::isfinite(pressure) || pressure <= 0.0
            || !std::isfinite(phaseEnthalpy)) {
            throw std::runtime_error(
                "Eulerian primary state is invalid during primitive recovery.");
        }

        for (int component = 0; component < 3; ++component) {
            const double momentum =
                phase.primary.momentum[(size_t)component]
                    .values()[(size_t)cell];
            if (!std::isfinite(momentum)) {
                throw std::runtime_error(
                    "Eulerian primary momentum is non-finite.");
            }
            const double velocity = momentum / mass;
            phase.primitive.velocity[(size_t)component]
                .values()[(size_t)cell] = velocity;
        }

        const double enthalpy = phaseEnthalpy / mass;
        const double temperature =
            Multiphase::phaseTemperatureFromEnthalpy(
                properties, enthalpy);
        if (!std::isfinite(enthalpy) || enthalpy <= 0.0
            || !std::isfinite(temperature) || temperature <= 0.0) {
            throw std::runtime_error(
                "Eulerian primary enthalpy produced invalid temperature.");
        }

        const double density = Multiphase::phaseDensity(
            properties, temperature, pressure);

        const double alpha = mass / density;
        if (!std::isfinite(density) || density <= 0.0
            || !std::isfinite(alpha) || alpha < 0.0 || alpha > 1.0) {
            throw std::runtime_error(
                "Eulerian thermo closure produced invalid primitive cache.");
        }
        phase.primitive.density.values()[(size_t)cell] = density;
        phase.primitive.alpha.values()[(size_t)cell] = alpha;
        phase.primitive.temperature.values()[(size_t)cell] = temperature;
        phase.primitive.enthalpy.values()[(size_t)cell] = enthalpy;
    }
}

void EulerianEulerianPhaseSystem::reconstructReferencePhaseMass() {
    if (phases_.empty()) {
        throw std::runtime_error(
            "Cannot close reference phase before setup.");
    }
    PhaseState& reference = phases_[referencePhase_];
    const auto& properties = phaseProperties(referencePhase_);
    for (int cell = 0; cell < geometry().TotalSize(); ++cell) {
        double independentAlpha = 0.0;
        for (size_t phase = 0; phase < phases_.size(); ++phase) {
            if (phase == referencePhase_) continue;
            independentAlpha += phases_[phase].primitive.alpha
                .values()[(size_t)cell];
        }
        const double alpha = 1.0 - independentAlpha;
        const double density =
            reference.primitive.density.values()[(size_t)cell];
        if (!std::isfinite(alpha) || alpha <= 0.0 || alpha > 1.0
            || !std::isfinite(density) || density <= 0.0) {
            throw std::runtime_error(
                "N-1 phaseMass closure produced invalid reference alpha.");
        }

        const double mass = alpha * density;
        for (int component = 0; component < 3; ++component) {
            const double velocity =
                reference.primitive.velocity[(size_t)component]
                    .values()[(size_t)cell];
            reference.primary.momentum[(size_t)component]
                .values()[(size_t)cell] = mass * velocity;
        }
        const double temperature =
            reference.primitive.temperature.values()[(size_t)cell];
        reference.primary.phaseMass.values()[(size_t)cell] = mass;
        reference.primary.phaseEnthalpy.values()[(size_t)cell] =
            mass * Multiphase::phaseSpecificEnthalpy(
                properties, temperature);
        reference.primitive.alpha.values()[(size_t)cell] = alpha;
    }
}

void EulerianEulerianPhaseSystem::recoverPrimitiveState() {
    for (size_t phase = 0; phase < phases_.size(); ++phase) {
        recoverPhase(phase);
    }
    reconstructReferencePhaseMass();
}

void EulerianEulerianPhaseSystem::computeInterphase(
        double dt,
        const std::vector<PhaseVectorField>& previousVelocity) {
    if (previousVelocity.size() != phases_.size()) {
        throw std::runtime_error(
            "Interphase previous-velocity workspace size is invalid.");
    }
    std::vector<Multiphase::PhaseProperties> properties;
    properties.reserve(properties_.size());
    for (const auto* phase : properties_) properties.push_back(*phase);
    Interphase::compute(
        {geometry(), config_, properties, phases_, previousVelocity, dt},
        sources_);
}

void EulerianEulerianPhaseSystem::validateState(
        const std::string& stage) const {
    for (int cell = 0; cell < geometry().TotalSize(); ++cell) {
        const double pressure = pressure_.values()[(size_t)cell];
        if (!std::isfinite(pressure) || pressure <= 0.0) {
            throw std::runtime_error(
                stage + ": shared pressure is invalid.");
        }
        double alphaSum = 0.0;
        for (size_t phaseIndex = 0;
             phaseIndex < phases_.size(); ++phaseIndex) {
            const PhaseState& phase = phases_[phaseIndex];
            const double mass =
                phase.primary.phaseMass.values()[(size_t)cell];
            const double alpha =
                phase.primitive.alpha.values()[(size_t)cell];
            const double density =
                phase.primitive.density.values()[(size_t)cell];
            const double temperature =
                phase.primitive.temperature.values()[(size_t)cell];
            if (!std::isfinite(mass) || mass <= 0.0
                || !std::isfinite(alpha) || alpha <= 0.0 || alpha > 1.0
                || !std::isfinite(density) || density <= 0.0
                || !std::isfinite(temperature) || temperature <= 0.0) {
                throw std::runtime_error(
                    stage + ": Eulerian phase state is invalid.");
            }
            const double massScale = std::max(std::abs(mass), 1.0);
            if (std::abs(mass - alpha * density)
                > 1.0e-11 * massScale) {
                throw std::runtime_error(
                    stage + ": phaseMass and primitive cache are inconsistent.");
            }
            for (int component = 0; component < 3; ++component) {
                const double velocity =
                    phase.primitive.velocity[(size_t)component]
                        .values()[(size_t)cell];
                const double momentum =
                    phase.primary.momentum[(size_t)component]
                        .values()[(size_t)cell];
                if (!std::isfinite(velocity)
                    || std::abs(momentum - mass * velocity)
                       > 1.0e-11 * std::max(std::abs(momentum), 1.0)) {
                    throw std::runtime_error(
                        stage + ": phase momentum cache is inconsistent.");
                }
            }
            const auto& properties = phaseProperties(phaseIndex);
            const double expectedEnthalpy = mass
                * Multiphase::phaseSpecificEnthalpy(
                    properties, temperature);
            if (std::abs(
                    phase.primary.phaseEnthalpy.values()[(size_t)cell]
                    - expectedEnthalpy)
                > 1.0e-10 * std::max(std::abs(expectedEnthalpy), 1.0)) {
                throw std::runtime_error(
                    stage + ": phase enthalpy cache is inconsistent.");
            }
            alphaSum += alpha;
        }
        if (std::abs(alphaSum - 1.0) > 1.0e-12) {
            throw std::runtime_error(
                stage + ": sum(alpha_k) differs from one.");
        }
    }
}

} // namespace SF::Physics::PhaseSystems
