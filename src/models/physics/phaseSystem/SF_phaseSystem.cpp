/// @file SF_phaseSystem.cpp
/// @brief PhaseSystem 状态、边界、源项与双欧拉物理实现。

#include "SF_phaseSystem.h"

#include "SF_eulerianEulerian.h"

#include <algorithm>
#include <stdexcept>

namespace SF::Physics::PhaseSystems {

void PhaseEquationSources::setupLike(
        const Field& field,
        const std::vector<PhaseState>& phases,
        const std::vector<std::array<size_t, 2>>& phasePairs) {
    mass.resize(phases.size());
    momentum.resize(phases.size());
    energy.resize(phases.size());
    wallBoilingMass.setupLike(field, "wallBoilingMass");
    wallBoilingConvectiveHeat.setupLike(
        field, "wallBoilingConvectiveHeat");
    wallBoilingQuenchingHeat.setupLike(
        field, "wallBoilingQuenchingHeat");
    wallBoilingEvaporativeHeat.setupLike(
        field, "wallBoilingEvaporativeHeat");
    wallBoilingDepartureDiameter.setupLike(
        field, "wallBoilingDepartureDiameter");
    interphaseMechanicalHeating.setupLike(
        field, "interphaseMechanicalHeating");
    transferLedger=InterphaseTransfer::TransferLedger(
        field.TotalSize(),(int)phases.size());
    for (size_t k = 0; k < phases.size(); ++k) {
        mass[k].setupLike(field, "interphaseMass." + phases[k].name);
        energy[k].setupLike(field, "interphaseEnergy." + phases[k].name);
        for (int d = 0; d < 3; ++d) {
            momentum[k][(size_t)d].setupLike(
                field, "interphaseMomentum." + phases[k].name);
        }
    }
    momentumCouplings.resize(phasePairs.size());
    for (size_t pair = 0; pair < phasePairs.size(); ++pair) {
        auto& coupling = momentumCouplings[pair];
        coupling.first = phasePairs[pair][0];
        coupling.second = phasePairs[pair][1];
        const std::string suffix = phases[coupling.first].name + "."
                                 + phases[coupling.second].name;
        coupling.coefficient.setupLike(
            field, "interphaseMomentumCoefficient." + suffix);
        coupling.dragCoefficient.setupLike(
            field, "interphaseDragCoefficient." + suffix);
    }
}

void PhaseEquationSources::clear() {
    transferLedger=InterphaseTransfer::TransferLedger(
        mass.front().TotalSize(),(int)mass.size());
    for (auto& f : mass) std::fill(f.values().begin(), f.values().end(), 0.0);
    for (auto& p : momentum)
        for (auto& f : p) std::fill(f.values().begin(), f.values().end(), 0.0);
    for (auto& f : energy) std::fill(f.values().begin(), f.values().end(), 0.0);
    std::fill(interphaseMechanicalHeating.values().begin(),
              interphaseMechanicalHeating.values().end(), 0.0);
    for (auto& coupling : momentumCouplings) {
        std::fill(coupling.coefficient.values().begin(),
                  coupling.coefficient.values().end(), 0.0);
        std::fill(coupling.dragCoefficient.values().begin(),
                  coupling.dragCoefficient.values().end(), 0.0);
    }
}

void PhaseEquationSources::clearExternalLedger() {
    for (ScalarField* field : {
             &wallBoilingMass, &wallBoilingConvectiveHeat,
             &wallBoilingQuenchingHeat, &wallBoilingEvaporativeHeat}) {
        std::fill(field->values().begin(), field->values().end(), 0.0);
    }
}

std::unique_ptr<PhaseSystem> makePhaseSystem(
        const Multiphase::MultiPhaseConfig& config) {
    if (Multiphase::isEulerianEulerianType(config.type)) {
        return std::make_unique<EulerianEulerianPhaseSystem>(config);
    }
    throw std::runtime_error(
        "PhaseSystem factory only owns the eulerianEulerian multi-velocity path; "
        "homogeneousMixture and volumeOfFluid remain on their dedicated paths.");
}

} // namespace SF::Physics::PhaseSystems
