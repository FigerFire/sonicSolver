#pragma once

/// @file SF_equationComposition.h
/// @brief Case-level model identifiers and equation composition input.

#include <string>
#include <unordered_set>
#include <vector>

namespace SF {

/// @brief Open provider registry. Adding a model does not edit a central enum.
class CompositionProviderRegistry {
public:
    void registerEquationOfState(std::string id) { equationOfState_.insert(std::move(id)); }
    void registerCaloricThermo(std::string id) { caloricThermo_.insert(std::move(id)); }
    void registerTransport(std::string id) { transport_.insert(std::move(id)); }
    void registerAlgorithm(std::string id) { algorithms_.insert(std::move(id)); }
    [[nodiscard]] bool hasEquationOfState(const std::string& id) const { return equationOfState_.count(id) != 0; }
    [[nodiscard]] bool hasCaloricThermo(const std::string& id) const { return caloricThermo_.count(id) != 0; }
    [[nodiscard]] bool hasTransport(const std::string& id) const { return transport_.count(id) != 0; }
    [[nodiscard]] bool hasAlgorithm(const std::string& id) const { return algorithms_.count(id) != 0; }
private:
    std::unordered_set<std::string> equationOfState_, caloricThermo_, transport_, algorithms_;
};

inline CompositionProviderRegistry builtinCompositionProviders() {
    CompositionProviderRegistry result;
    for (const auto* id : {"perfectGas", "rhoConst"}) result.registerEquationOfState(id);
    for (const auto* id : {"hConst", "janaf"}) result.registerCaloricThermo(id);
    for (const auto* id : {"const", "sutherland"}) result.registerTransport(id);
    for (const auto* id : {"Explicit", "SIMPLE", "PISO", "PIMPLE"}) result.registerAlgorithm(id);
    return result;
}

struct ThermoDynamicsSelection {
    std::string equationOfState;
    std::string thermo;
    std::string transport;
    double constantDensity = 0.0;
};

/// @brief Independently named mathematical-system input object.
struct EquationSystemInstanceConfig {
    std::string name;
    std::string type;
};

struct EquationCompositionConfig {
    bool declared = false;
    std::vector<std::string> equations;
    ThermoDynamicsSelection thermoDynamics;
    std::string algorithm;
    int outerCorrectors = 1;
    int pressureCorrectors = 1;
    int nonOrthogonalCorrectors = 0;
    std::vector<EquationSystemInstanceConfig> instances;
};

} // namespace SF
