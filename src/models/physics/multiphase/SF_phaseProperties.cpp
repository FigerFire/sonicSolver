/// @file SF_phaseProperties.cpp
/// @brief 多相配置、物性与模型一致性校验实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_phaseProperties.h"
#include "SF_phaseChange.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <set>
#include <stdexcept>
#include <utility>

namespace SF {
namespace Physics {
namespace Multiphase {

namespace {

std::string compactLower(std::string s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        if (c == ' ' || c == '_' || c == '-') continue;
        out.push_back((char)std::tolower((unsigned char)c));
    }
    return out;
}

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

bool phaseNeedsPressureBased(const PhaseProperties& phase) {
    const std::string thermo = normalizeModelType(phase.thermoModel);
    return thermo == "perfectliquid" || thermo == "solid";
}

PhaseRole resolvedRole(const MultiPhaseConfig& config,
                       const std::string& phaseName) {
    if (const PhaseProperties* phase = findPhase(config, phaseName)) {
        if (phase->role != PhaseRole::Unknown) return phase->role;
        const PhaseRole thermoRole = inferRoleFromThermoModel(
            phase->thermoModel);
        if (thermoRole != PhaseRole::Unknown) return thermoRole;
    }
    return inferRoleFromName(phaseName);
}

} // namespace

std::string normalizeModelType(std::string type) {
    return compactLower(std::move(type));
}

std::string normalizePhaseName(std::string name) {
    return compactLower(std::move(name));
}

PhaseRole parsePhaseRole(std::string role) {
    const std::string r = compactLower(std::move(role));
    if (r == "liquid" || r == "water" || r == "positive") {
        return PhaseRole::Liquid;
    }
    if (r == "gas" || r == "air" || r == "vapor" || r == "vapour"
        || r == "negative") {
        return PhaseRole::Gas;
    }
    return PhaseRole::Unknown;
}

std::string phaseRoleName(PhaseRole role) {
    switch (role) {
    case PhaseRole::Liquid: return "liquid";
    case PhaseRole::Gas: return "gas";
    default: return "unknown";
    }
}

bool isLevelSetType(const std::string& type) {
    return normalizeModelType(type) == "levelset";
}

bool isMixtureType(const std::string& type) {
    const std::string t = normalizeModelType(type);
    return t == "mixture"
        || t == "singlevelocitymixture"
        || t == "singletemperaturesinglevelocitymixture"
        || t == "singlevelocitysingletemperaturemixture"
        || t == "mixturelike"
        || t == "alpha";
}

bool isHomogeneousType(const std::string& type) {
    const std::string t = normalizeModelType(type);
    return t == "homogeneous" || t == "homogeneousmultiphase"
        || t == "singlepressuresinglevelocitysingletemperature";
}

bool isEulerianEulerianType(const std::string& type) {
    const std::string t = normalizeModelType(type);
    return t == "eulerianeulerian" || t == "twomomentum";
}

bool isThermalType(const std::string& type) {
    const std::string t = normalizeModelType(type);
    return t == "thermal"
        || t == "singletemperature"
        || t == "heatconduction"
        || t == "temperature";
}

MultiPhaseSolverKind parseMultiPhaseSolverKind(std::string solver) {
    const std::string raw = solver;
    const std::string s = normalizeModelType(std::move(solver));
    if (s.empty() || s == "auto") return MultiPhaseSolverKind::Auto;
    if (s == "densitybased" || s == "density") {
        return MultiPhaseSolverKind::DensityBased;
    }
    if (s == "pressurebased" || s == "pressure") {
        return MultiPhaseSolverKind::PressureBased;
    }
    throw std::runtime_error(
        "multiPhase.solver supports auto, densityBased, or pressureBased; got '"
        + raw + "'.");
}

std::string multiPhaseSolverKindName(MultiPhaseSolverKind solver) {
    switch (solver) {
    case MultiPhaseSolverKind::DensityBased: return "densityBased";
    case MultiPhaseSolverKind::PressureBased: return "pressureBased";
    default: return "auto";
    }
}

MultiPhaseSolverKind resolvedMultiPhaseSolverKind(
        const MultiPhaseConfig& config) {
    if (config.solver != MultiPhaseSolverKind::Auto) return config.solver;
    const std::string phaseChangeModel =
        normalizeModelType(config.phaseChange.model);
    if (config.phaseChange.enabled
        && (phaseChangeModel == "rpi"
            || phaseChangeModel == "wallboiling")) {
        return MultiPhaseSolverKind::PressureBased;
    }
    for (const PhaseProperties& phase : config.phases) {
        if (phaseNeedsPressureBased(phase)) {
            return MultiPhaseSolverKind::PressureBased;
        }
    }
    return MultiPhaseSolverKind::DensityBased;
}

std::vector<std::string> densityBasedConservedVariables(
        const MultiPhaseConfig& config) {
    std::vector<std::string> variables = {
        "rho", "rhoU", "rhoV", "rhoW", "rhoE"
    };
    for (const PhaseProperties& phase : config.phases) {
        if (normalizeModelType(phase.thermoModel) == "perfectgas") {
            variables.push_back("rhoY." + phase.name);
        }
    }
    for (const std::string& alpha : config.mixture.phaseFractionFields) {
        const std::string variable = "rhoAlpha." + alpha;
        if (std::find(variables.begin(), variables.end(), variable)
            == variables.end()) {
            variables.push_back(variable);
        }
    }
    return variables;
}

const PhaseProperties* findPhase(const MultiPhaseConfig& config,
                                 const std::string& name) {
    const std::string wanted = normalizePhaseName(name);
    for (const PhaseProperties& phase : config.phases) {
        if (normalizePhaseName(phase.name) == wanted) return &phase;
    }
    return nullptr;
}

double phaseSign(const MultiPhaseConfig& config,
                 const std::string& phaseName) {
    const PhaseRole role = resolvedRole(config, phaseName);
    if (role == PhaseRole::Liquid) return 1.0;
    if (role == PhaseRole::Gas) return -1.0;
    throw std::runtime_error(
        "multiPhase: phase '" + phaseName
        + "' has unknown role. Declare [phase." + phaseName
        + "] role = \"liquid\" or role = \"gas\".");
}

} // namespace Multiphase
} // namespace Physics
} // namespace SF
