/// @file SF_phaseState.cpp
/// @brief PhaseSystem 状态、边界、源项与双欧拉物理实现。

#include "SF_phaseState.h"

#include <cmath>
#include <stdexcept>

namespace SF::Physics::PhaseSystems {

void PhasePrimaryState::setupLike(const Field& field,
                                  const std::string& phaseName) {
    phaseMass.setupLike(field, "phaseMass." + phaseName);
    phaseEnthalpy.setupLike(field, "phaseEnthalpy." + phaseName);
    for (int d = 0; d < 3; ++d) {
        const std::string component(1, "xyz"[d]);
        momentum[(size_t)d].setupLike(
            field, "momentum" + component + "." + phaseName);
    }
}

void PhasePrimitiveCache::setupLike(const Field& field,
                                    const std::string& phaseName) {
    alpha.setupLike(field, "alpha." + phaseName);
    density.setupLike(field, "rho." + phaseName);
    temperature.setupLike(field, "T." + phaseName);
    enthalpy.setupLike(field, "h." + phaseName);
    for (int d = 0; d < 3; ++d) {
        const std::string component(1, "xyz"[d]);
        velocity[(size_t)d].setupLike(
            field, "U" + component + "." + phaseName);
    }
}

void PhaseState::setupLike(const Field& field,
                           const std::string& phaseName) {
    name = phaseName;
    primary.setupLike(field, phaseName);
    primitive.setupLike(field, phaseName);
}

void PhaseState::commitPrimitiveToPrimary(
        const Multiphase::PhaseProperties& properties) {
    for (int n = 0; n < primitive.alpha.TotalSize(); ++n) {
        const double a = primitive.alpha.values()[(size_t)n];
        const double rho = primitive.density.values()[(size_t)n];
        const double T = primitive.temperature.values()[(size_t)n];
        if (!std::isfinite(a) || !std::isfinite(rho) || rho <= 0.0
            || !std::isfinite(T) || T <= 0.0) {
            throw std::runtime_error(
                "PhaseState cannot commit invalid primitive state.");
        }
        const double mass = a * rho;
        primary.phaseMass.values()[(size_t)n] = mass;
        for (int d = 0; d < 3; ++d) {
            const double u = primitive.velocity[(size_t)d].values()[(size_t)n];
            if (!std::isfinite(u)) {
                throw std::runtime_error("PhaseState found non-finite velocity.");
            }
            primary.momentum[(size_t)d].values()[(size_t)n] = mass * u;
        }
        primary.phaseEnthalpy.values()[(size_t)n] =
            mass * Multiphase::phaseSpecificEnthalpy(properties, T);
    }
}

} // namespace SF::Physics::PhaseSystems
