/// @file SF_phaseValidation.cpp
/// @brief 多相配置、物性与模型一致性校验实现。

#include "SF_phaseProperties.h"
#include "SF_phaseChange.h"
#include "SF_phaseValidation.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <set>
#include <stdexcept>

namespace SF {
namespace Physics {
namespace Multiphase {

namespace {

PhaseRole inferRoleFromName(const std::string& name) {
    const std::string n = normalizePhaseName(name);
    if (n == "water" || n == "liquid" || n == "fluid" || n == "oil") {
        return PhaseRole::Liquid;
    }
    if (n == "air" || n == "gas" || n == "vapor" || n == "vapour"
        || n == "steam" || n == "bubble") {
        return PhaseRole::Gas;
    }
    return PhaseRole::Unknown;
}

PhaseRole inferRoleFromThermoModel(const std::string& model) {
    const std::string m = normalizeModelType(model);
    if (m == "perfectliquid" || m == "liquid" || m == "incompressibleliquid") {
        return PhaseRole::Liquid;
    }
    if (m == "perfectgas" || m == "idealgas" || m == "gas") {
        return PhaseRole::Gas;
    }
    return PhaseRole::Unknown;
}

} // namespace

void validateMultiPhaseConfig(const MultiPhaseConfig& config,
                              const std::string& context) {
    if (!config.enabled) return;

    const bool levelSetType = isLevelSetType(config.type);
    const bool mixtureType = isMixtureType(config.type);
    const bool thermalType = isThermalType(config.type);
    const bool homogeneousType = isHomogeneousType(config.type);
    const bool eulerianEulerianType = isEulerianEulerianType(config.type);
    const MultiPhaseSolverKind solverKind =
        resolvedMultiPhaseSolverKind(config);

    if (!levelSetType && !mixtureType && !thermalType && !homogeneousType
        && !eulerianEulerianType) {
        throw std::runtime_error(
            context + ": [multiPhase].type supports \"level set\" or "
            "\"mixture\", \"eulerianEulerian\", and thermophysicalProperties "
            "supports \"thermal\"; got '"
            + config.type + "'.");
    }

    if (homogeneousType) {
        if (!std::isfinite(config.initialPressure) || config.initialPressure <= 0.0) {
            throw std::runtime_error(context + ": homogeneous FluidStateModel requires finite initialPressure > 0.");
        }
        if (!config.temperature.enabled || !std::isfinite(config.temperature.defaultValue)
            || config.temperature.defaultValue <= 0.0) {
            throw std::runtime_error(context + ": homogeneous FluidStateModel requires [temperature] enabled with a positive defaultValue.");
        }
        double alphaSum = 0.0;
        for (const PhaseProperties& phase : config.phases) {
            if (normalizeModelType(phase.thermoModel) != "stiffenedgas"
                || !std::isfinite(phase.volumeFraction) || phase.volumeFraction < 0.0) {
                throw std::runtime_error(context + ": homogeneous baseline requires stiffenedGas phases with explicit non-negative volumeFraction.");
            }
            alphaSum += phase.volumeFraction;
            if (!phase.species.empty()) {
                double ySum = 0.0;
                for (const SpeciesProperties& species : phase.species) {
                    if (species.name.empty() || !std::isfinite(species.massFraction)
                        || species.massFraction < 0.0) {
                        throw std::runtime_error(context + ": homogeneous species require name and non-negative massFraction.");
                    }
                    ySum += species.massFraction;
                }
                if (std::abs(ySum - 1.0) > 1.0e-12) {
                    throw std::runtime_error(context + ": each phase's species massFraction values must sum exactly to one.");
                }
            }
        }
        if (config.phases.size() < 2 || std::abs(alphaSum - 1.0) > 1.0e-12) {
            throw std::runtime_error(context + ": homogeneous baseline requires at least two phases whose volumeFraction values sum exactly to one.");
        }
    }

    if (eulerianEulerianType) {
        if (config.eulerianEulerian.phaseNames.size() < 2) {
            throw std::runtime_error(
                context + ": eulerianEulerian requires at least two phases.");
        }
        if (findPhase(config, config.eulerianEulerian.referencePhase) == nullptr) {
            throw std::runtime_error(
                context + ": eulerianEulerian referencePhase is not declared.");
        }
        if (config.eulerianEulerian.phasePairs.empty()) {
            throw std::runtime_error(
                context + ": interphaseModels.phasePairs must register at "
                "least one phase pair.");
        }
        if (config.eulerianEulerian.radialCoordinate < 0
            || config.eulerianEulerian.radialCoordinate > 2) {
            throw std::runtime_error(
                context + ": Eulerian axisymmetric radialCoordinate "
                "must be x, y or z.");
        }
        std::set<std::string> pairNames;
        std::set<std::string> pairKeys;
        for (const auto& pair : config.eulerianEulerian.phasePairs) {
            if (pair.name.empty()
                || !pairNames.insert(normalizePhaseName(pair.name)).second) {
                throw std::runtime_error(
                    context + ": Eulerian phase-pair names must be non-empty "
                    "and unique.");
            }
            if (findPhase(config, pair.continuousPhase) == nullptr
                || findPhase(config, pair.dispersedPhase) == nullptr
                || normalizePhaseName(pair.continuousPhase)
                    == normalizePhaseName(pair.dispersedPhase)) {
                throw std::runtime_error(
                    context + ": phase pair '" + pair.name
                    + "' needs two distinct declared phases.");
            }
            std::array<std::string, 2> names{
                normalizePhaseName(pair.continuousPhase),
                normalizePhaseName(pair.dispersedPhase)};
            std::sort(names.begin(), names.end());
            if (!pairKeys.insert(names[0] + "|" + names[1]).second) {
                throw std::runtime_error(
                    context + ": duplicate unordered phase pair '"
                    + pair.name + "'.");
            }
            try {
                (void)phasePairDiameter(config, pair);
            } catch (const std::exception& error) {
                throw std::runtime_error(
                    context + ": " + error.what());
            }
            if (!std::isfinite(pair.liftCoefficient)
                || !std::isfinite(pair.virtualMassCoefficient)
                || pair.virtualMassCoefficient < 0.0) {
                throw std::runtime_error(
                    context + ": phase pair '" + pair.name
                    + "' has invalid drag/lift/virtual-mass coefficients.");
            }
            const std::string wallModel =
                normalizeModelType(pair.wallLubricationModel);
            if (wallModel != "none" && wallModel != "antal") {
                throw std::runtime_error(
                    context + ": phase pair '" + pair.name
                    + "' wallLubrication supports none or Antal.");
            }
            if (wallModel == "antal"
                && (pair.wallPatch.empty()
                    || !std::isfinite(pair.wallLubricationC1)
                    || !std::isfinite(pair.wallLubricationC2))) {
                throw std::runtime_error(
                    context + ": Antal wallLubrication for pair '"
                    + pair.name + "' needs wallPatch and finite C1/C2.");
            }
            const std::string dispersionModel =
                normalizeModelType(pair.turbulentDispersionModel);
            if (dispersionModel != "none"
                && dispersionModel != "gdb"
                && dispersionModel != "favreaveraged") {
                throw std::runtime_error(
                    context + ": phase pair '" + pair.name
                    + "' turbulentDispersion supports none or GDB.");
            }
            if (dispersionModel != "none"
                && (!std::isfinite(pair.turbulentDispersionCoefficient)
                    || pair.turbulentDispersionCoefficient < 0.0
                    || !std::isfinite(pair.turbulentKinematicViscosity)
                    || pair.turbulentKinematicViscosity <= 0.0
                    || !std::isfinite(pair.turbulentSchmidtNumber)
                    || pair.turbulentSchmidtNumber <= 0.0)) {
                throw std::runtime_error(
                    context + ": turbulentDispersion for pair '" + pair.name
                    + "' needs C0>=0, nuT>0 and sigmaAlpha>0.");
            }
            const std::string tensionModel =
                normalizeModelType(pair.surfaceTensionModel);
            if (tensionModel != "none"
                && (tensionModel != "constant"
                    || !std::isfinite(pair.surfaceTension)
                    || pair.surfaceTension <= 0.0)) {
                throw std::runtime_error(
                    context + ": surfaceTension for pair '" + pair.name
                    + "' supports none or constant sigma>0.");
            }
        }
        for (const PhaseProperties& phase : config.phases) {
            if (!std::isfinite(phase.density) || phase.density <= 0.0
                || !std::isfinite(phase.volumeFraction)
                || phase.volumeFraction < 0.0 || phase.volumeFraction > 1.0
                || !std::isfinite(phase.Cv) || phase.Cv <= 0.0
                || !std::isfinite(phase.viscosity) || phase.viscosity < 0.0
                || !std::isfinite(phase.thermalConductivity)
                || phase.thermalConductivity < 0.0) {
                throw std::runtime_error(
                    context + ": Eulerian-Eulerian phase '" + phase.name
                    + "' needs SI rho/alpha/Cv/mu/k values.");
            }
            if (!std::isfinite(phase.residualAlpha)
                || phase.residualAlpha < 0.0
                || phase.residualAlpha >= 1.0) {
                throw std::runtime_error(
                    context + ": phase '" + phase.name
                    + "' residualAlpha must be in [0,1).");
            }
        }
        SF::Physics::PhaseChange::validateConfig(config, context);
    }

    std::set<std::string> phaseNames;
    for (const PhaseProperties& phase : config.phases) {
        if (phase.name.empty()) {
            throw std::runtime_error(context + ": phase name cannot be empty.");
        }
        const std::string normalizedName = normalizePhaseName(phase.name);
        if (!phaseNames.insert(normalizedName).second) {
            throw std::runtime_error(
                context + ": duplicate phase declaration '" + phase.name + "'.");
        }
        if (phase.role == PhaseRole::Unknown
            && inferRoleFromThermoModel(phase.thermoModel)
                == PhaseRole::Unknown
            && inferRoleFromName(phase.name) == PhaseRole::Unknown) {
            throw std::runtime_error(
                context + ": phase '" + phase.name
                + "' needs model = \"perfectLiquid\" or "
                "model = \"perfectGas\".");
        }
        auto validateLaw = [&](const TemperaturePropertyLaw& law,
                               double fallback,
                               const std::string& property) {
            const std::string model = normalizeModelType(law.model);
            if (model.empty() || model == "constant") return;
            if (model != "polynomialtemperature"
                || law.coefficients.empty()
                || !std::isfinite(law.referenceTemperature)
                || !std::isfinite(law.minimumTemperature)
                || !std::isfinite(law.maximumTemperature)
                || law.minimumTemperature <= 0.0
                || law.maximumTemperature <= law.minimumTemperature) {
                throw std::runtime_error(
                    context + ": phase '" + phase.name + "' " + property
                    + " polynomialTemperature needs finite Tref, "
                      "validTemperature (Tmin Tmax), and coefficients.");
            }
            for (double temperature : {
                     law.minimumTemperature,
                     0.5 * (law.minimumTemperature
                            + law.maximumTemperature),
                     law.maximumTemperature}) {
                (void)temperaturePropertyValue(
                    law, fallback, temperature,
                    property + "(" + phase.name + ")");
            }
        };
        validateLaw(
            phase.densityLaw, phase.density, "density");
        validateLaw(
            phase.viscosityLaw, phase.viscosity, "dynamicViscosity");
        validateLaw(
            phase.specificHeatLaw, phase.specificHeat, "specificHeat");
        validateLaw(
            phase.thermalConductivityLaw,
            phase.thermalConductivity, "thermalConductivity");
    }

    auto validateTemperature = [&]() {
        if (!config.temperature.enabled) {
            if (thermalType) {
                throw std::runtime_error(
                    context + ": thermal model requires [temperature].enabled true.");
            }
            return;
        }
        if (!std::isfinite(config.temperature.diffusivity)
            || config.temperature.diffusivity < 0.0) {
            throw std::runtime_error(
                context + ": [temperature].diffusivity must be finite and >= 0.");
        }
        if (std::isfinite(config.temperature.defaultValue)
            && config.temperature.defaultValue <= 0.0) {
            throw std::runtime_error(
                context + ": [temperature].defaultValue must be > 0 when set.");
        }
        if (!std::isfinite(config.temperature.lowerBound)
            || !std::isfinite(config.temperature.upperBound)
            || config.temperature.lowerBound > config.temperature.upperBound) {
            throw std::runtime_error(
                context + ": [temperature].bounds must be finite and ordered.");
        }
        for (const auto& condition : config.temperature.initialConditions) {
            if (condition.name.empty()) {
                throw std::runtime_error(
                    context + ": IC [T] set name cannot be empty.");
            }
            if (condition.type != FIXED_VALUE) {
                throw std::runtime_error(
                    context + ": IC [T] set '" + condition.name
                    + "' supports only FIXED_VALUE.");
            }
            if (!std::isfinite(condition.value) || condition.value <= 0.0) {
                throw std::runtime_error(
                    context + ": IC [T] set '" + condition.name
                    + "' value must be finite and > 0.");
            }
        }
        for (const auto& condition : config.temperature.boundaryConditions) {
            if (condition.name.empty()) {
                throw std::runtime_error(
                    context + ": BC [T] set name cannot be empty.");
            }
            if (condition.type == FIXED_VALUE
                && (!std::isfinite(condition.value)
                    || condition.value <= 0.0)) {
                throw std::runtime_error(
                    context + ": BC [T] set '" + condition.name
                    + "' value must be finite and > 0.");
            }
        }
    };

    validateLevelSetConfig(config, context);

    if (mixtureType) {
        if (config.phases.size() < 2) {
            throw std::runtime_error(
                context + ": mixture model requires at least two explicit phases.");
        }
        for (const auto& phase : config.phases) {
            const std::string thermo = normalizeModelType(phase.thermoModel);
            if (thermo != "perfectliquid" && thermo != "perfectgas"
                && thermo != "stiffenedgas" && thermo != "solid") {
                throw std::runtime_error(
                    context + ": mixture phase '" + phase.name
                    + "' must declare an explicit supported EOS model.");
            }
            const PhaseRole modelRole = inferRoleFromThermoModel(
                phase.thermoModel);
            const PhaseRole explicitRole =
                phase.role == PhaseRole::Unknown
                    ? inferRoleFromName(phase.name) : phase.role;
            if (modelRole != PhaseRole::Unknown
                && explicitRole != PhaseRole::Unknown
                && explicitRole != modelRole) {
                throw std::runtime_error(
                    context + ": phase '" + phase.name
                    + "' role conflicts with thermo model '"
                    + phase.thermoModel + "'.");
            }
        }
        if (!std::isfinite(config.alpha.defaultValue)) {
            throw std::runtime_error(
                context + ": [alpha].defaultValue must be finite.");
        }
        if (!std::isfinite(config.alpha.diffusivity)
            || config.alpha.diffusivity < 0.0) {
            throw std::runtime_error(
                context + ": [alpha].diffusivity must be finite and >= 0.");
        }
        if (!std::isfinite(config.alpha.lowerBound)
            || !std::isfinite(config.alpha.upperBound)
            || config.alpha.lowerBound > config.alpha.upperBound) {
            throw std::runtime_error(
                context + ": [alpha].bounds must be finite and ordered.");
        }
        for (const auto& condition : config.alpha.initialConditions) {
            if (condition.name.empty()) {
                throw std::runtime_error(
                    context + ": IC [alpha] set name cannot be empty.");
            }
            if (condition.type != FIXED_VALUE) {
                throw std::runtime_error(
                    context + ": IC [alpha] set '" + condition.name
                    + "' supports only FIXED_VALUE.");
            }
            if (!std::isfinite(condition.value)) {
                throw std::runtime_error(
                    context + ": IC [alpha] set '" + condition.name
                    + "' value must be finite.");
            }
        }
        for (const auto& condition : config.alpha.boundaryConditions) {
            if (condition.name.empty()) {
                throw std::runtime_error(
                    context + ": BC [alpha] set name cannot be empty.");
            }
            if (condition.type == FIXED_VALUE
                && !std::isfinite(condition.value)) {
                throw std::runtime_error(
                    context + ": BC [alpha] set '" + condition.name
                    + "' value must be finite.");
            }
        }

        validateTemperature();

        SF::Physics::PhaseChange::validateConfig(config, context);
    }

    if (thermalType) {
        validateTemperature();
        if (config.phaseChange.enabled
            && normalizeModelType(config.phaseChange.model) != "none") {
            throw std::runtime_error(
                context + ": standalone thermal model does not support "
                "[phaseChange]; use a mixture model for phase change.");
        }
    }
}

} // namespace Multiphase
} // namespace Physics
} // namespace SF
