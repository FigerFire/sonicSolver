/// @file SF_homogeneousPhaseChange.cpp
/// @brief 守恒变量 FluidStateModel 与工厂实现。

#include "SF_homogeneousPhaseChange.h"

#include "methods/numerics/structured/SF_structured.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace SF::Physics::FluidStateModel {
namespace {

const Multiphase::PhaseProperties& phaseByRole(
        const Multiphase::MultiPhaseConfig& config,
        Multiphase::PhaseRole role) {
    for (const auto& phase : config.phases) {
        if (phase.role == role) return phase;
    }
    throw std::runtime_error("Homogeneous phase change requires explicit liquid and gas phase roles.");
}

int transferComponent(
        const HomogeneousMultiphaseStateModel& equations,
        const Multiphase::PhaseProperties& phase,
        const Multiphase::PhaseChangeOptions& phaseChange,
        const char* selectionKey) {
    const std::string species = PhaseChange::selection(phaseChange, selectionKey, "");
    if (!species.empty()) {
        const int index = equations.componentIndex(phase.name, species);
        if (index < 0) {
            throw std::runtime_error(
                "Phase-change transfer species '" + species
                + "' is not declared in phase '" + phase.name + "'.");
        }
        return index;
    }
    return equations.uniqueComponentIndex(phase.name);
}

} // namespace

HomogeneousPhaseChange::HomogeneousPhaseChange(
        Multiphase::MultiPhaseConfig config,
        std::shared_ptr<const HomogeneousMultiphaseStateModel> equations)
    : config_(std::move(config)), equations_(std::move(equations)) {
    if (!equations_) throw std::runtime_error("HomogeneousPhaseChange requires FluidStateModel.");
}

void HomogeneousPhaseChange::initialize(const Field& field) {
    if (!config_.phaseChange.enabled) return;
    const auto& liquid = phaseByRole(config_, Multiphase::PhaseRole::Liquid);
    const auto& vapor = phaseByRole(config_, Multiphase::PhaseRole::Gas);
    liquidComponent_ = transferComponent(*equations_, liquid,
                                         config_.phaseChange, "liquidSpecies");
    vaporComponent_ = transferComponent(*equations_, vapor,
                                        config_.phaseChange, "vaporSpecies");
    liquidPhase_ = equations_->findPhaseIndex(liquid.name);
    vaporPhase_ = equations_->findPhaseIndex(vapor.name);
    if (liquidPhase_ < 0 || vaporPhase_ < 0 || liquidPhase_ == vaporPhase_) {
        throw std::runtime_error(
            "Homogeneous phase change could not map distinct liquid/vapor phases.");
    }
    trackedPhase_ = equations_->findPhaseIndex(config_.alpha.phaseName);
    if (trackedPhase_ < 0) {
        throw std::runtime_error("alpha phase '" + config_.alpha.phaseName
                                 + "' is not in Homogeneous FluidStateModel.");
    }
    alpha_.setupLike(field, "alpha." + config_.alpha.phaseName, 0.0);
    temperature_.setupLike(field, "T", 0.0);
    refreshDerived(field);
}

void HomogeneousPhaseChange::refreshDerived(const Field& field) {
    for (int k = 0; k < field.MZ(); ++k)
        for (int j = 0; j < field.MY(); ++j)
            for (int i = 0; i < field.MX(); ++i) {
                const auto state = field.thermodynamicState(i, j, k);
                if (trackedPhase_ >= static_cast<int>(state.volumeFraction.size())
                    || !std::isfinite(state.temperature)
                    || state.temperature <= 0.0) {
                    throw std::runtime_error("Homogeneous phase-change derived state is invalid.");
                }
                alpha_(i,j,k) = state.volumeFraction[(size_t)trackedPhase_];
                temperature_(i,j,k) = state.temperature;
            }
}

void HomogeneousPhaseChange::assemble(Field& field, Residual& residual, double dt) {
    if (!config_.phaseChange.enabled) return;
    refreshDerived(field);
    PhaseChange::PhaseChangeRateResult rates =
        PhaseChange::computeRates(field, alpha_, temperature_, config_, dt);
    ledger_ = InterphaseTransfer::TransferLedger(
        field.TotalSize(), equations_->phaseCount());
    rates.diagnostics.minMdot = 1.0e300;
    rates.diagnostics.maxMdot = -1.0e300;
    Math::forFluidInterior(field, [&](int i, int j, int k) {
        const int id = alpha_.getIdx(i,j,k);
        double mdot = rates.mdot[(size_t)id];
        if (!std::isfinite(mdot)) {
            throw std::runtime_error("Homogeneous phase-change mdot is non-finite.");
        }
        const auto state = field.thermodynamicState(i,j,k);
        if (liquidPhase_ >= static_cast<int>(state.phaseMass.size())
            || vaporPhase_ >= static_cast<int>(state.phaseMass.size())) {
            throw std::runtime_error(
                "Homogeneous phase-change closure did not return phase masses.");
        }
        const double liquidLimit =
            state.phaseMass[(size_t)liquidPhase_] / dt;
        const double vaporLimit =
            state.phaseMass[(size_t)vaporPhase_] / dt;
        if (mdot > liquidLimit || -mdot > vaporLimit) {
            throw std::runtime_error(
                "Homogeneous phase-change source exhausts its donor phase at cell ("
                + std::to_string(i) + "," + std::to_string(j) + ","
                + std::to_string(k) + "): mdot="
                + std::to_string(mdot) + " kg/(m3 s), liquidLimit="
                + std::to_string(liquidLimit) + ", vaporLimit="
                + std::to_string(vaporLimit)
                + ". Reduce maxDeltaT or the phase-change coefficient.");
        }
        rates.mdot[(size_t)id] = mdot;
        rates.diagnostics.minMdot = std::min(rates.diagnostics.minMdot, mdot);
        rates.diagnostics.maxMdot = std::max(rates.diagnostics.maxMdot, mdot);
        if (mdot > 0.0) {
            ledger_.addInternalMassTransfer(
                id, liquidPhase_, vaporPhase_, mdot);
        } else if (mdot < 0.0) {
            ledger_.addInternalMassTransfer(
                id, vaporPhase_, liquidPhase_, -mdot);
        }
    });
    ledger_.validate();

    std::vector<int> routes((size_t)equations_->phaseCount(), -1);
    for (int phase = 0; phase < equations_->phaseCount(); ++phase) {
        routes[(size_t)phase] = equations_->firstComponentIndex(phase);
    }
    routes[(size_t)liquidPhase_] = liquidComponent_;
    routes[(size_t)vaporPhase_] = vaporComponent_;
    InterphaseTransfer::SourceBundle sources(
        field.TotalSize(), equations_->variableCount());
    equations_->mapTransfersToSources(ledger_, routes, sources);
    Math::forFluidInterior(field, [&](int i, int j, int k) {
        const int id = field.getIdx(i,j,k);
        for (int variable = 0; variable < field.NVar(); ++variable) {
            residual.source(i,j,k,variable) += sources(id, variable);
        }
    });
    diagnostics_ = rates.diagnostics;
}

} // namespace SF::Physics::FluidStateModel
