#pragma once

/// @file SF_equationComposition.h
/// @brief Case-level model identifiers and equation composition input.

#include <string>
#include <memory>
#include <map>
#include <limits>
#include <cmath>
#include <stdexcept>
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

/// @brief 归一化后冻结的物理闭合；kernel 投影不是另一个可编辑物性来源。
struct ThermophysicalContract {
    ThermoDynamicsSelection selection;
    double gamma = std::numeric_limits<double>::quiet_NaN();
    double gasConstant = std::numeric_limits<double>::quiet_NaN();
    double dynamicViscosity = 0.0;
    double prandtl = 0.72;
    bool native = false;
    std::string source;
    std::map<std::string,std::string> provenance;
    [[nodiscard]] double cv() const { return gasConstant/(gamma-1.0); }
    [[nodiscard]] double cp() const { return gamma*gasConstant/(gamma-1.0); }
    [[nodiscard]] double conductivity() const { return dynamicViscosity*cp()/prandtl; }
    void validate() const {
        if (!std::isfinite(dynamicViscosity)||dynamicViscosity<0.0||!std::isfinite(prandtl)||prandtl<=0.0)
            throw std::runtime_error(source+": const transport requires finite mu>=0 and Pr>0.");

        if (selection.equationOfState=="perfectGas") {
            if (!std::isfinite(gamma)||gamma<=1.0||!std::isfinite(gasConstant)||gasConstant<=0.0)
                throw std::runtime_error(source+": perfectGas requires finite gamma>1 and R>0.");
            if (!std::isfinite(cp())||!std::isfinite(cv())||!std::isfinite(conductivity()))
                throw std::runtime_error(source+": derived cp/cv/k must be finite.");
        } else if (selection.equationOfState=="rhoConst") {
            if (!std::isfinite(selection.constantDensity)||selection.constantDensity<=0.0)
                throw std::runtime_error(source+": rhoConst requires finite rho>0.");
        } else throw std::runtime_error(source+": unknown equationOfState '"+selection.equationOfState+"'.");
        if (!selection.thermo.empty()&&selection.thermo!="hConst")
            throw std::runtime_error(source+": Unsupported caloric provider '"+selection.thermo+"'; implemented: hConst.");
        if (!selection.transport.empty()&&selection.transport!="const")
            throw std::runtime_error(source+": Unsupported transport provider '"+selection.transport+"'; implemented: const.");
    }
};

/// @brief Independently named mathematical-system input object.
struct EquationSystemInstanceConfig {
    std::string name;
    std::string type;
};

struct EquationCompositionConfig {
    bool declared = false;
    /// Solution membership is independent of mathematical role and dependency activation.
    bool stateDeclared = false;
    std::vector<std::string> solutionVariables;
    std::string stateSelectionOrigin;
    /// Compatibility adapters lower old equation-source requests once, before composition.
    bool compatibilityPressureConstraint = false;
    std::vector<std::string> equations;
    ThermoDynamicsSelection thermoDynamics;
    std::string algorithm;
    int outerCorrectors = 1;
    int pressureCorrectors = 1;
    int nonOrthogonalCorrectors = 0;
    std::vector<EquationSystemInstanceConfig> instances;
};

} // namespace SF
