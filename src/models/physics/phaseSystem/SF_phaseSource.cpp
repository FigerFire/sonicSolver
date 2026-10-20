/// @file SF_phaseSource.cpp
/// @brief PhaseSystem 状态、边界、源项与双欧拉物理实现。

#include "SF_phaseSource.h"

#include <cmath>
#include <stdexcept>

namespace SF::Physics::PhaseSystems {

void PhaseSourceRegistry::add(
        std::unique_ptr<PhaseEquationSource> source) {
    if (!source) {
        throw std::runtime_error(
            "PhaseSourceRegistry cannot register a null provider.");
    }
    sources_.push_back(std::move(source));
}

void PhaseSourceRegistry::assemble(
        PhaseSystem& system,
        double dt) const {
    if (!std::isfinite(dt) || dt <= 0.0) {
        throw std::runtime_error(
            "PhaseSourceRegistry requires finite positive dt.");
    }
    system.sources().clearExternalLedger();
    for (const auto& source : sources_) {
        source->add(system, dt, system.sources());
    }
    system.sources().transferLedger.validate();
}

std::vector<std::string> PhaseSourceRegistry::names() const {
    std::vector<std::string> result;
    result.reserve(sources_.size());
    for (const auto& source : sources_) result.push_back(source->name());
    return result;
}

} // namespace SF::Physics::PhaseSystems
