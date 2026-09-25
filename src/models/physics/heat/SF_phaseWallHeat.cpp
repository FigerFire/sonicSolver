/// @file SF_phaseWallHeat.cpp
/// @brief 壁面热源与相能量耦合实现。

#include "SF_phaseWallHeat.h"

#include "SF_closure.h"
#include "SF_heatFluxPartition.h"
#include "SF_phaseChange.h"
#include "SF_wallMapping.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace SF::Physics::PhaseSystems {
namespace {

size_t uniqueLiquidPhase(const PhaseSystem& system) {
    size_t result = system.phases().size();
    for (size_t phase = 0; phase < system.phases().size(); ++phase) {
        if (system.phaseProperties(phase).role
            != Multiphase::PhaseRole::Liquid) {
            continue;
        }
        if (result != system.phases().size()) {
            throw std::runtime_error(
                "WallHeat requires one unambiguous liquid target phase.");
        }
        result = phase;
    }
    if (result == system.phases().size()) {
        throw std::runtime_error(
            "WallHeat requires one phase with role liquid.");
    }
    return result;
}

size_t uniqueGasPhase(const PhaseSystem& system) {
    size_t result = system.phases().size();
    for (size_t phase = 0; phase < system.phases().size(); ++phase) {
        if (system.phaseProperties(phase).role
            != Multiphase::PhaseRole::Gas) {
            continue;
        }
        if (result != system.phases().size()) {
            throw std::runtime_error(
                "WallBoiling requires one unambiguous gas phase.");
        }
        result = phase;
    }
    if (result == system.phases().size()) {
        throw std::runtime_error(
            "WallBoiling requires one phase with role gas.");
    }
    return result;
}

size_t configuredPhase(
        const PhaseSystem& system,
        Multiphase::PhaseRole role) {
    const auto& names = system.config().phaseChange.phaseNames;
    if (names.empty()) {
        return role == Multiphase::PhaseRole::Liquid
            ? uniqueLiquidPhase(system) : uniqueGasPhase(system);
    }
    size_t result = system.phases().size();
    for (const std::string& name : names) {
        for (size_t phase = 0; phase < system.phases().size(); ++phase) {
            if (Multiphase::normalizePhaseName(
                    system.phases()[phase].name)
                    != Multiphase::normalizePhaseName(name)
                || system.phaseProperties(phase).role != role) {
                continue;
            }
            if (result != system.phases().size()) {
                throw std::runtime_error(
                    "phaseChange phases contain an ambiguous phase role.");
            }
            result = phase;
        }
    }
    if (result == system.phases().size()) {
        throw std::runtime_error(
            "phaseChange phases do not contain the required "
            + Multiphase::phaseRoleName(role) + " phase.");
    }
    return result;
}

class PhaseWallHeatSource final : public PhaseEquationSource {
public:
    explicit PhaseWallHeatSource(std::vector<WallHeatSetting> settings)
        : settings_(std::move(settings)) {
        if (settings_.empty()) {
            throw std::runtime_error(
                "WallHeat is enabled but no wall heat setting was loaded.");
        }
        for (const auto& setting : settings_) {
            if (setting.patch.empty() || setting.patch == "all"
                || !std::isfinite(setting.heatFlux)
                || !setting.coupleEnergy) {
                throw std::runtime_error(
                    "WallHeat needs a named patch, finite heatFlux and "
                    "coupleEnergy=true.");
            }
            if (setting.rangeCoordinate >= 0
                && (!std::isfinite(setting.rangeMinimum)
                    || !std::isfinite(setting.rangeMaximum)
                    || setting.rangeMaximum <= setting.rangeMinimum)) {
                throw std::runtime_error(
                    "WallHeat coordinate range must be finite and ordered.");
            }
            if (std::isfinite(setting.wallTemperature)
                && setting.wallTemperature <= 0.0) {
                throw std::runtime_error(
                    "WallHeat wallTemperature must be positive.");
            }
        }
    }

    std::string name() const override { return "WallHeat"; }

    void add(
            const PhaseSystem& system,
            double dt,
            PhaseEquationSources& sources) const override {
        const size_t liquid = configuredPhase(
            system, Multiphase::PhaseRole::Liquid);
        const size_t vapor = configuredPhase(
            system, Multiphase::PhaseRole::Gas);
        const auto& liquidState = system.phases()[liquid];
        const bool rpi = PhaseChange::normalizeModel(
            system.config().phaseChange.model) == "rpi";
        if (cachedGeometry_ != &system.geometry()
            || cachedSamples_.size() != settings_.size()) {
            cachedGeometry_ = &system.geometry();
            cachedSamples_.clear();
            cachedSamples_.reserve(settings_.size());
            for (const auto& setting : settings_) {
                Multiphase::PhaseProperties liquidProperties =
                    system.phaseProperties(liquid);
                Multiphase::PhaseProperties vaporProperties =
                    system.phaseProperties(vapor);
                const PhaseChange::ModelContext mappingContext{
                    system.geometry(), liquidState.primitive.alpha,
                    liquidState.primitive.temperature, system.config(),
                    liquidProperties, vaporProperties, liquidProperties,
                    true, dt};
                cachedSamples_.push_back(
                    PhaseChange::RPI::mapWallPatch(
                        mappingContext, setting));
            }
        }
        for (size_t settingIndex = 0;
             settingIndex < settings_.size(); ++settingIndex) {
            const auto& setting = settings_[settingIndex];
            Multiphase::PhaseProperties liquidProperties =
                system.phaseProperties(liquid);
            Multiphase::PhaseProperties vaporProperties =
                system.phaseProperties(vapor);
            const auto& samples = cachedSamples_[settingIndex];
            for (const auto& sample : samples) {
                const int cell = system.geometry().getIdx(
                    sample.fluidI, sample.fluidJ, sample.fluidK);
                if (!rpi) {
                    const double source =
                        setting.heatFlux * sample.areaOverVolume;
                    if (!std::isfinite(source)) {
                        throw std::runtime_error(
                            "WallHeat produced non-finite volumetric source.");
                    }
                    sources.energy[liquid].values()[(size_t)cell] += source;
                    continue;
                }

                liquidProperties.density = liquidState.primitive.density
                    .values()[(size_t)cell];
                vaporProperties.density = system.phases()[vapor]
                    .primitive.density.values()[(size_t)cell];
                const double localLiquidTemperature =
                    liquidState.primitive.temperature.values()[(size_t)cell];
                const double vaporTemperature =
                    system.phases()[vapor].primitive.temperature
                        .values()[(size_t)cell];
                liquidProperties.viscosity =
                    Multiphase::phaseViscosity(
                        liquidProperties, localLiquidTemperature);
                liquidProperties.specificHeat =
                    Multiphase::phaseSpecificHeat(
                        liquidProperties, localLiquidTemperature);
                liquidProperties.thermalConductivity =
                    Multiphase::phaseThermalConductivity(
                        liquidProperties, localLiquidTemperature);
                vaporProperties.viscosity =
                    Multiphase::phaseViscosity(
                        vaporProperties, vaporTemperature);
                vaporProperties.specificHeat =
                    Multiphase::phaseSpecificHeat(
                        vaporProperties, vaporTemperature);
                vaporProperties.thermalConductivity =
                    Multiphase::phaseThermalConductivity(
                        vaporProperties, vaporTemperature);
                double liquidSpeedSquared = 0.0;
                for (int component = 0; component < 3; ++component) {
                    const double velocity = liquidState.primitive.velocity
                        [(size_t)component].values()[(size_t)cell];
                    liquidSpeedSquared += velocity * velocity;
                }
                const double localSpeed = std::sqrt(liquidSpeedSquared);
                const double localFrictionVelocity = std::sqrt(
                    liquidProperties.viscosity * localSpeed
                    / (liquidProperties.density
                       * sample.wallNormalDistance));
                PhaseChange::ModelContext localContext{
                    system.geometry(), liquidState.primitive.alpha,
                    liquidState.primitive.temperature, system.config(),
                    liquidProperties, vaporProperties, liquidProperties,
                    true, dt,
                    system.sharedPressure().values()[(size_t)cell],
                    localSpeed, sample.wallNormalDistance,
                    localFrictionVelocity};
                const double wallTemperature =
                    std::isfinite(setting.wallTemperature)
                    ? setting.wallTemperature
                    : liquidState.primitive.temperature(
                          sample.wallI, sample.wallJ, sample.wallK);
                const double liquidTemperature =
                    liquidState.primitive.temperature(
                        sample.fluidI, sample.fluidJ, sample.fluidK);
                const auto balance =
                    PhaseChange::RPI::evaluateWallHeatBalance(
                        localContext, liquidTemperature,
                        wallTemperature, setting.heatFlux);
                const auto& closure = balance.closure;
                const auto& heat = balance.heat;
                const double maximumVaporFraction =
                    PhaseChange::coefficient(
                        system.config().phaseChange,
                        "maximumWallVaporFraction", -1.0);
                const double vaporFraction = system.phases()[vapor]
                    .primitive.alpha.values()[(size_t)cell];
                if (!std::isfinite(maximumVaporFraction)
                    || maximumVaporFraction <= 0.0
                    || maximumVaporFraction >= 1.0) {
                    throw std::runtime_error(
                        "RPI wall ledger requires explicit "
                        "maximumWallVaporFraction in (0,1).");
                }
                if (vaporFraction >= maximumVaporFraction) {
                    throw std::runtime_error(
                        "RPI wall vapor fraction reached the configured "
                        "nucleate-boiling validity limit; dryout/CHF "
                        "transition is not implemented.");
                }
                const double tolerance = PhaseChange::coefficient(
                    system.config().phaseChange,
                    "heatFluxBalanceTolerance", -1.0);
                if (!std::isfinite(tolerance) || tolerance <= 0.0) {
                    throw std::runtime_error(
                        "RPI wall ledger requires explicit positive "
                        "heatFluxBalanceTolerance.");
                }
                const double imbalance =
                    std::abs(heat.total - setting.heatFlux);
                if (imbalance > tolerance * setting.heatFlux) {
                    std::ostringstream message;
                    message << "RPI wall ledger heat imbalance on patch '"
                            << setting.patch << "': qWall="
                            << setting.heatFlux << ", qConv="
                            << heat.convective << ", qQuench="
                            << heat.quenching << ", qEvap="
                            << heat.evaporative << ", tolerance="
                            << tolerance << ".";
                    throw std::runtime_error(message.str());
                }
                const double mdot =
                    closure.wallMassFlux * sample.areaOverVolume;
                if (!std::isfinite(mdot) || mdot <= 0.0) {
                    throw std::runtime_error(
                        "RPI wall ledger produced invalid volumetric mass rate.");
                }
                const double hLiquid = liquidState.primitive.enthalpy
                    .values()[(size_t)cell];
                const double qConv = heat.convective
                    * sample.areaOverVolume;
                const double qQuench = heat.quenching
                    * sample.areaOverVolume;
                const double qEvap = heat.evaporative
                    * sample.areaOverVolume;
                sources.transferLedger.addInternalMassTransfer(
                    cell,(int)liquid,(int)vapor,mdot);
                std::array<double,3> carriedMomentum{};
                sources.mass[liquid].values()[(size_t)cell] -= mdot;
                sources.mass[vapor].values()[(size_t)cell] += mdot;
                for (int component = 0; component < 3; ++component) {
                    const double momentum = mdot
                        * liquidState.primitive.velocity[(size_t)component]
                            .values()[(size_t)cell];
                    carriedMomentum[(size_t)component]=momentum;
                    sources.momentum[liquid][(size_t)component]
                        .values()[(size_t)cell] -= momentum;
                    sources.momentum[vapor][(size_t)component]
                        .values()[(size_t)cell] += momentum;
                }
                sources.transferLedger.addInternalMomentumTransfer(
                    cell,(int)liquid,(int)vapor,carriedMomentum);
                sources.transferLedger.addInternalEnergyTransfer(
                    cell,(int)liquid,(int)vapor,mdot*hLiquid);
                sources.transferLedger.addWallEnergyInput(
                    cell,(int)liquid,qConv+qQuench);
                sources.transferLedger.addWallEnergyInput(
                    cell,(int)vapor,qEvap);
                sources.energy[liquid].values()[(size_t)cell] +=
                    qConv + qQuench - mdot * hLiquid;
                sources.energy[vapor].values()[(size_t)cell] +=
                    mdot * hLiquid + qEvap;
                sources.wallBoilingMass.values()[(size_t)cell] += mdot;
                sources.wallBoilingConvectiveHeat.values()[(size_t)cell]
                    += qConv;
                sources.wallBoilingQuenchingHeat.values()[(size_t)cell]
                    += qQuench;
                sources.wallBoilingEvaporativeHeat.values()[(size_t)cell]
                    += qEvap;
                sources.wallBoilingDepartureDiameter.values()[(size_t)cell]
                    = closure.departureDiameter;
            }
        }
    }

private:
    std::vector<WallHeatSetting> settings_;
    mutable const Field* cachedGeometry_ = nullptr;
    mutable std::vector<std::vector<PhaseChange::RPI::WallCellSample>>
        cachedSamples_;
};

} // namespace

std::unique_ptr<PhaseEquationSource> makePhaseWallHeatSource(
        std::vector<WallHeatSetting> settings) {
    return std::make_unique<PhaseWallHeatSource>(std::move(settings));
}

} // namespace SF::Physics::PhaseSystems
