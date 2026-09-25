/// @file SF_temperatureProperties.cpp
/// @brief 多相配置、物性与模型一致性校验实现。

#include "SF_phaseProperties.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace SF::Physics::Multiphase {
namespace {

bool polynomial(const TemperaturePropertyLaw& law) {
    return normalizeModelType(law.model) == "polynomialtemperature";
}

void requireTemperatureRange(const TemperaturePropertyLaw& law,
                             double temperature,
                             const std::string& propertyName) {
    if (!std::isfinite(temperature)
        || temperature < law.minimumTemperature
        || temperature > law.maximumTemperature) {
        throw std::runtime_error(
            propertyName + " temperature " + std::to_string(temperature)
            + " K is outside the declared polynomial range ["
            + std::to_string(law.minimumTemperature) + ", "
            + std::to_string(law.maximumTemperature) + "] K.");
    }
}

double polynomialValue(const TemperaturePropertyLaw& law,
                       double temperature) {
    double value = 0.0;
    const double theta = temperature - law.referenceTemperature;
    for (auto it = law.coefficients.rbegin();
         it != law.coefficients.rend(); ++it) {
        value = value * theta + *it;
    }
    return value;
}

double polynomialIntegral(const TemperaturePropertyLaw& law,
                          double lower,
                          double upper) {
    const double lo = lower - law.referenceTemperature;
    const double hi = upper - law.referenceTemperature;
    double integral = 0.0;
    double loPower = lo;
    double hiPower = hi;
    for (size_t order = 0; order < law.coefficients.size(); ++order) {
        integral += law.coefficients[order]
            * (hiPower - loPower) / (double)(order + 1);
        loPower *= lo;
        hiPower *= hi;
    }
    return integral;
}

double constantCp(const PhaseProperties& phase) {
    const double cp = phase.specificHeat > 0.0
        ? phase.specificHeat : phase.Cv + phase.gasConstant;
    if (!std::isfinite(cp) || cp <= 0.0) {
        throw std::runtime_error(
            "Phase '" + phase.name + "' requires finite positive Cp.");
    }
    return cp;
}

} // namespace

double temperaturePropertyValue(const TemperaturePropertyLaw& law,
                                double fallback,
                                double temperature,
                                const std::string& propertyName) {
    const std::string model = normalizeModelType(law.model);
    if (model.empty() || model == "constant") {
        if (!std::isfinite(fallback)) {
            throw std::runtime_error(propertyName + " is non-finite.");
        }
        return fallback;
    }
    if (!polynomial(law)) {
        throw std::runtime_error(
            propertyName + " supports constant or polynomialTemperature.");
    }
    requireTemperatureRange(law, temperature, propertyName);
    if (law.coefficients.empty()) {
        throw std::runtime_error(
            propertyName + " polynomial has no coefficients.");
    }
    const double value = polynomialValue(law, temperature);
    if (!std::isfinite(value) || value <= 0.0) {
        throw std::runtime_error(
            propertyName + " polynomial produced a non-positive value.");
    }
    return value;
}

double phaseDensity(const PhaseProperties& phase,
                    double temperature,
                    double pressure) {
    if (polynomial(phase.densityLaw)) {
        return temperaturePropertyValue(
            phase.densityLaw, phase.density, temperature,
            "rho(" + phase.name + ")");
    }
    const std::string model = normalizeModelType(phase.thermoModel);
    if (model == "perfectgas") {
        if (!std::isfinite(pressure) || pressure <= 0.0
            || !std::isfinite(phase.gasConstant)
            || phase.gasConstant <= 0.0) {
            throw std::runtime_error(
                "perfectGas density needs positive p and gasConstant.");
        }
        return pressure / (phase.gasConstant * temperature);
    }
    if (model == "perfectliquid" || model == "constant") {
        return temperaturePropertyValue(
            phase.densityLaw, phase.density, temperature,
            "rho(" + phase.name + ")");
    }
    throw std::runtime_error(
        "Unsupported thermoModel '" + phase.thermoModel + "'.");
}

double phaseViscosity(const PhaseProperties& phase, double temperature) {
    return temperaturePropertyValue(
        phase.viscosityLaw, phase.viscosity, temperature,
        "mu(" + phase.name + ")");
}

double phaseSpecificHeat(const PhaseProperties& phase, double temperature) {
    return temperaturePropertyValue(
        phase.specificHeatLaw, constantCp(phase), temperature,
        "Cp(" + phase.name + ")");
}

double phaseThermalConductivity(const PhaseProperties& phase,
                                double temperature) {
    return temperaturePropertyValue(
        phase.thermalConductivityLaw, phase.thermalConductivity, temperature,
        "k(" + phase.name + ")");
}

double phaseSpecificEnthalpy(const PhaseProperties& phase,
                             double temperature) {
    if (!polynomial(phase.specificHeatLaw)) {
        return constantCp(phase) * temperature;
    }
    const auto& law = phase.specificHeatLaw;
    requireTemperatureRange(law, temperature, "h(" + phase.name + ")");
    const double cpReference = polynomialValue(
        law, law.referenceTemperature);
    const double referenceEnthalpy =
        cpReference * law.referenceTemperature;
    const double enthalpy = referenceEnthalpy + polynomialIntegral(
        law, law.referenceTemperature, temperature);
    if (!std::isfinite(enthalpy) || enthalpy <= 0.0) {
        throw std::runtime_error(
            "Cp polynomial produced invalid phase enthalpy.");
    }
    return enthalpy;
}

double phaseTemperatureFromEnthalpy(const PhaseProperties& phase,
                                    double enthalpy) {
    if (!std::isfinite(enthalpy) || enthalpy <= 0.0) {
        throw std::runtime_error("Phase enthalpy must be finite and positive.");
    }
    if (!polynomial(phase.specificHeatLaw)) {
        return enthalpy / constantCp(phase);
    }
    const auto& law = phase.specificHeatLaw;
    double lower = law.minimumTemperature;
    double upper = law.maximumTemperature;
    const double hLower = phaseSpecificEnthalpy(phase, lower);
    const double hUpper = phaseSpecificEnthalpy(phase, upper);
    if (enthalpy < hLower || enthalpy > hUpper) {
        throw std::runtime_error(
            "Phase '" + phase.name + "' enthalpy "
            + std::to_string(enthalpy)
            + " J/kg is outside Cp(T) integral range ["
            + std::to_string(hLower) + ", "
            + std::to_string(hUpper) + "] J/kg.");
    }
    for (int iteration = 0; iteration < 80; ++iteration) {
        const double middle = 0.5 * (lower + upper);
        const double hMiddle = phaseSpecificEnthalpy(phase, middle);
        if (hMiddle < enthalpy) lower = middle;
        else upper = middle;
    }
    return 0.5 * (lower + upper);
}

double phasePairDiameter(const MultiPhaseConfig& config,
                         const PhasePairModelOptions& pair) {
    if (std::isfinite(pair.particleDiameter)
        && pair.particleDiameter > 0.0) {
        return pair.particleDiameter;
    }
    const PhaseProperties* dispersed = findPhase(
        config, pair.dispersedPhase);
    if (!dispersed
        || normalizeModelType(dispersed->diameterModel) != "constant"
        || !std::isfinite(dispersed->diameter)
        || dispersed->diameter <= 0.0) {
        throw std::runtime_error(
            "Phase pair '" + pair.name
            + "' needs a positive pair diameter or constant dispersed-phase "
              "diameter.");
    }
    return dispersed->diameter;
}

double phasePairSurfaceTension(const MultiPhaseConfig& config,
                               const std::string& first,
                               const std::string& second) {
    const std::string a = normalizePhaseName(first);
    const std::string b = normalizePhaseName(second);
    for (const auto& pair : config.eulerianEulerian.phasePairs) {
        const std::string continuous =
            normalizePhaseName(pair.continuousPhase);
        const std::string dispersed =
            normalizePhaseName(pair.dispersedPhase);
        if (!((continuous == a && dispersed == b)
              || (continuous == b && dispersed == a))) {
            continue;
        }
        if (normalizeModelType(pair.surfaceTensionModel) != "constant"
            || !std::isfinite(pair.surfaceTension)
            || pair.surfaceTension <= 0.0) {
            throw std::runtime_error(
                "Phase pair '" + pair.name
                + "' requires constant positive surface tension.");
        }
        return pair.surfaceTension;
    }
    throw std::runtime_error(
        "No surface-tension model is registered for phase pair '"
        + first + " and " + second + "'.");
}

} // namespace SF::Physics::Multiphase
