/// @file SF_nativeDecode.cpp
/// @brief Native typed description decoding.

#include "SF_compatibility.h"
#include "app/application/model/SF_configParser.h"
#include "core/interfaces/SF_log.h"
#include "SF_cellFaceMesh.h"
#include "core/mesh/SF_dimension.h"
#include "SF_meshGen.h"
#include "SF_phaseChange.h"
#include "SF_multiphase.h"
#include "SF_config.h"
#include "private/SF_casePath.h"
#include "private/SF_valueParser.h"
#include "private/SF_sourceParserUtils.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>
#include <sys/stat.h>

namespace SF::IOPrivate {

using P = Model::Parameters;

std::string token(const P& value) {
    if (value.is_string()) {
        const std::string s = value.get<std::string>();
        if (s.find('/') != std::string::npos
            || s.find(' ') != std::string::npos) {
            return "\"" + s + "\"";
        }
        return s;
    }
    if (value.is_array()) {
        std::string s = "(";
        for (const auto& x : value) s += token(x) + " ";
        return s + ")";
    }
    return value.dump();
}

namespace {

bool modelEnabled(const P* value,std::initializer_list<const char*> wrappers={}) {
    if (!value || !value->is_object()) return false;
    std::optional<bool> enabled;
    const auto read=[&](const P& block) {
        if(block.contains("enable"))throw std::runtime_error("Optional model uses reserved 'enabled', not 'enable'.");
        if(!block.contains("enabled"))return;
        if(!block.at("enabled").is_boolean())throw std::runtime_error("Model enabled must be boolean.");
        const bool selected=block.at("enabled").get<bool>();
        if(enabled && *enabled!=selected)throw std::runtime_error("Conflicting explicit model enabled at outer and nested blocks.");
        enabled=selected;
    };
    read(*value);
    for(const auto* name:wrappers)if(value->contains(name) && value->at(name).is_object())read(value->at(name));
    return enabled.value_or(true);
}

void collectKeyMatches(const P* node, const std::string& key,
                       std::vector<std::pair<const P*, int>>& out) {
    if (node == nullptr || !node->is_object()) return;
    for (auto it = node->begin(); it != node->end(); ++it) {
        const bool isObject = it.value().is_object();
        if (it.key() == key) {
            out.emplace_back(&it.value(), isObject ? 2 : 1);
        }
        if (isObject) collectKeyMatches(&it.value(), key, out);
    }
}

std::string jsonValueAfterKey(const P* node, const std::string& key) {
    std::vector<std::pair<const P*, int>> matches;
    collectKeyMatches(node, key, matches);
    if (matches.empty() || matches.front().second != 1) return "";
    return token(*matches.front().first);
}

const P* jsonBlockAfterKey(const P* node, const std::string& key) {
    std::vector<std::pair<const P*, int>> matches;
    collectKeyMatches(node, key, matches);
    for (const auto& m : matches) {
        if (m.second == 2) return m.first;
    }
    return nullptr;
}

std::vector<std::string> jsonValuesAfterKey(const P* node,
                                            const std::string& key) {
    std::vector<std::pair<const P*, int>> matches;
    collectKeyMatches(node, key, matches);
    std::vector<std::string> out;
    out.reserve(matches.size());
    for (const auto& m : matches) {
        if (m.second == 1) out.push_back(token(*m.first));
    }
    return out;
}

std::vector<std::pair<std::string, const P*>> jsonChildBlocks(const P* node) {
    std::vector<std::pair<std::string, const P*>> out;
    if (node == nullptr || !node->is_object()) return out;
    for (auto it = node->begin(); it != node->end(); ++it) {
        if (it.value().is_object()) {
            out.emplace_back(it.key(), &it.value());
        }
    }
    return out;
}

} // namespace

} // namespace SF::IOPrivate

namespace SF::IOPrivate::PhasePropertiesReader {

using Config = Physics::Multiphase::MultiPhaseConfig;




void parseFieldModels(
        const P* text,
        const P* multiPhase,
        Config& config) {
auto parseLevelSet = [&](const P* block) {
    if (!block) return;
    std::string v = jsonValueAfterKey(block, "advectionOrder");
    if (!v.empty()) {
        config.levelSet.advectionOrder =
            valueIntValue(v, config.levelSet.advectionOrder);
    }
    v = jsonValueAfterKey(block, "advectionTimeScheme");
    if (v.empty()) v = jsonValueAfterKey(block, "timeScheme");
    if (!v.empty()) {
        fatalCaseConfig(
            "[levelSet].advectionTimeScheme/timeScheme has moved to "
            "solvers/numerics.yaml: time { default ...; }. The main flow "
            "and interface must use one explicit tableau.");
    }
    v = jsonValueAfterKey(block, "wenoEpsilon");
    if (!v.empty()) {
        config.levelSet.wenoEpsilon =
            valueDoubleValue(v, config.levelSet.wenoEpsilon);
    }
    v = jsonValueAfterKey(block, "wenoPower");
    if (!v.empty()) {
        config.levelSet.wenoPower =
            valueDoubleValue(v, config.levelSet.wenoPower);
    }
    v = jsonValueAfterKey(block, "advectFluidCellsOnly");
    if (!v.empty()) {
        const std::string key = valueLower(valueUnquote(v));
        if (key != "true" && key != "false" && key != "yes" && key != "no"
            && key != "on" && key != "off" && key != "1" && key != "0") {
            fatalCaseConfig(
                "[levelSet].advectFluidCellsOnly must be an explicit boolean.");
        }
        config.levelSet.advectFluidCellsOnly = valueBoolValue(v, false);
        config.levelSet.advectFluidCellsOnlyDeclared = true;
    }
    v = jsonValueAfterKey(block, "reinitializationOrder");
    if (!v.empty()) {
        config.levelSet.reinitializationOrder =
            valueIntValue(v, config.levelSet.reinitializationOrder);
    }
    v = jsonValueAfterKey(block, "reinitializationSteps");
    if (!v.empty()) {
        config.levelSet.reinitializationSteps =
            valueIntValue(v, config.levelSet.reinitializationSteps);
    }
    v = jsonValueAfterKey(block, "pseudoTimeStep");
    if (!v.empty()) {
        config.levelSet.pseudoTimeStep =
            valueDoubleValue(v, config.levelSet.pseudoTimeStep);
    }
    v = jsonValueAfterKey(block, "signSmoothingFactor");
    if (!v.empty()) {
        config.levelSet.signSmoothingFactor =
            valueDoubleValue(v, config.levelSet.signSmoothingFactor);
    }
    v = jsonValueAfterKey(block, "geometryGradientTolerance");
    if (!v.empty()) {
        config.levelSet.geometryGradientTolerance =
            valueDoubleValue(v, config.levelSet.geometryGradientTolerance);
    }
    v = jsonValueAfterKey(block, "interfaceThickness");
    if (!v.empty()) {
        config.levelSet.interfaceThickness =
            valueDoubleValue(v, config.levelSet.interfaceThickness);
    }
    v = jsonValueAfterKey(block, "narrowBandWidth");
    if (!v.empty()) {
        config.levelSet.narrowBandWidth =
            valueDoubleValue(v, config.levelSet.narrowBandWidth);
    }
    v = jsonValueAfterKey(block, "surfaceTension");
    if (!v.empty()) {
        config.levelSet.surfaceTension =
            valueDoubleValue(v, config.levelSet.surfaceTension);
    }
    v = jsonValueAfterKey(block, "surfaceTensionModel");
    if (!v.empty()) {
        config.levelSet.surfaceTensionModel = valueUnquote(v);
        config.levelSet.surfaceTensionModelDeclared = true;
    }
    v = jsonValueAfterKey(block, "massCorrection");
    if (!v.empty()) {
        const std::string key = valueLower(valueUnquote(v));
        if (key != "true" && key != "false" && key != "yes" && key != "no"
            && key != "on" && key != "off" && key != "1" && key != "0") {
            fatalCaseConfig(
                "[levelSet].massCorrection must be an explicit boolean.");
        }
        config.levelSet.massCorrection = valueBoolValue(v, false);
        config.levelSet.massCorrectionDeclared = true;
    }
};
parseLevelSet(multiPhase);
parseLevelSet(jsonBlockAfterKey(text, "levelSet"));

auto parseMixture = [&](const P* block) {
    if (!block) return;
    std::string v = jsonValueAfterKey(block, "phases");
    if (!v.empty()) config.mixture.phaseNames = valueWordList(v);
    for (const std::string& raw :
         jsonValuesAfterKey(block, "phaseFraction")) {
        for (const std::string& fieldName : valueWordList(raw)) {
            if (!fieldName.empty()
                && std::find(config.mixture.phaseFractionFields.begin(),
                             config.mixture.phaseFractionFields.end(),
                             fieldName)
                       == config.mixture.phaseFractionFields.end()) {
                config.mixture.phaseFractionFields.push_back(fieldName);
            }
        }
    }
    if (!config.mixture.phaseFractionFields.empty()) {
        config.alpha.fieldName = config.mixture.phaseFractionFields.front();
        const std::string first = config.alpha.fieldName;
        const size_t dot = first.find_last_of('.');
        if (dot != std::string::npos && dot + 1 < first.size()) {
            config.alpha.phaseName = first.substr(dot + 1);
        }
    }
};
parseMixture(multiPhase);
const P* mixtureBlock = jsonBlockAfterKey(text, "mixture");
if (!mixtureBlock) mixtureBlock = jsonBlockAfterKey(text, "mixtrue");
parseMixture(mixtureBlock);

auto parseAlpha = [&](const P* block) {
    if (!block) return;
    std::string v = jsonValueAfterKey(block, "fieldName");
    if (v.empty()) v = jsonValueAfterKey(block, "name");
    if (!v.empty()) config.alpha.fieldName = valueUnquote(v);
    v = jsonValueAfterKey(block, "phase");
    if (v.empty()) v = jsonValueAfterKey(block, "alphaPhase");
    if (!v.empty()) config.alpha.phaseName = valueUnquote(v);
    v = jsonValueAfterKey(block, "defaultAlpha");
    if (v.empty()) v = jsonValueAfterKey(block, "defaultValue");
    if (v.empty()) v = jsonValueAfterKey(block, "value");
    if (!v.empty()) {
        config.alpha.defaultValue =
            valueDoubleValue(v, config.alpha.defaultValue);
    }
    v = jsonValueAfterKey(block, "transport");
    if (v.empty()) v = jsonValueAfterKey(block, "transportEnabled");
    if (v.empty()) v = jsonValueAfterKey(block, "enabled");
    if (!v.empty()) {
        config.alpha.transportEnabled =
            valueBoolValue(v, config.alpha.transportEnabled);
    }
    v = jsonValueAfterKey(block, "diffusion");
    if (v.empty()) v = jsonValueAfterKey(block, "diffusionEnabled");
    if (!v.empty()) {
        config.alpha.diffusionEnabled =
            valueBoolValue(v, config.alpha.diffusionEnabled);
    }
    v = jsonValueAfterKey(block, "diffusivity");
    if (v.empty()) v = jsonValueAfterKey(block, "D");
    if (v.empty()) v = jsonValueAfterKey(block, "Dalpha");
    if (!v.empty()) {
        config.alpha.diffusivity =
            valueDoubleValue(v, config.alpha.diffusivity);
        config.alpha.diffusionEnabled = config.alpha.diffusivity > 0.0;
    }
    v = jsonValueAfterKey(block, "bounded");
    if (!v.empty()) {
        config.alpha.bounded =
            valueBoolValue(v, config.alpha.bounded);
    }
    v = jsonValueAfterKey(block, "boundMode");
    if (v.empty()) v = jsonValueAfterKey(block, "boundsMode");
    if (!v.empty()) config.alpha.boundMode = valueUnquote(v);
    v = jsonValueAfterKey(block, "lowerBound");
    if (v.empty()) v = jsonValueAfterKey(block, "min");
    if (!v.empty()) {
        config.alpha.lowerBound =
            valueDoubleValue(v, config.alpha.lowerBound);
    }
    v = jsonValueAfterKey(block, "upperBound");
    if (v.empty()) v = jsonValueAfterKey(block, "max");
    if (!v.empty()) {
        config.alpha.upperBound =
            valueDoubleValue(v, config.alpha.upperBound);
    }
    v = jsonValueAfterKey(block, "bounds");
    const std::vector<double> bounds = valueNumbers(v);
    if (bounds.size() >= 2) {
        config.alpha.lowerBound = bounds[0];
        config.alpha.upperBound = bounds[1];
    }
};
parseAlpha(multiPhase);
parseAlpha(jsonBlockAfterKey(text, "alpha"));

auto parseTemperature = [&](const P* block) {
    if (!block) return;
    std::string v = jsonValueAfterKey(block, "enabled");
    if (v.empty()) v = jsonValueAfterKey(block, "enable");
    if (v.empty()) v = jsonValueAfterKey(block, "solve");
    if (!v.empty()) {
        config.temperature.enabled =
            valueBoolValue(v, config.temperature.enabled);
    }
    v = jsonValueAfterKey(block, "fieldName");
    if (v.empty()) v = jsonValueAfterKey(block, "name");
    if (!v.empty()) config.temperature.fieldName = valueUnquote(v);
    v = jsonValueAfterKey(block, "defaultValue");
    if (v.empty()) v = jsonValueAfterKey(block, "defaultTemperature");
    if (v.empty()) v = jsonValueAfterKey(block, "T0");
    if (!v.empty()) {
        config.temperature.defaultValue =
            valueDoubleValue(v, config.temperature.defaultValue);
    }
    v = jsonValueAfterKey(block, "transport");
    if (v.empty()) v = jsonValueAfterKey(block, "transportEnabled");
    if (!v.empty()) {
        config.temperature.transportEnabled =
            valueBoolValue(v, config.temperature.transportEnabled);
    }
    v = jsonValueAfterKey(block, "diffusion");
    if (v.empty()) v = jsonValueAfterKey(block, "diffusionEnabled");
    if (!v.empty()) {
        config.temperature.diffusionEnabled =
            valueBoolValue(v, config.temperature.diffusionEnabled);
    }
    v = jsonValueAfterKey(block, "diffusivity");
    if (v.empty()) v = jsonValueAfterKey(block, "D");
    if (v.empty()) v = jsonValueAfterKey(block, "DT");
    if (!v.empty()) {
        config.temperature.diffusivity =
            valueDoubleValue(v, config.temperature.diffusivity);
        config.temperature.diffusionEnabled =
            config.temperature.diffusivity > 0.0;
    }
    v = jsonValueAfterKey(block, "bounded");
    if (!v.empty()) {
        config.temperature.bounded =
            valueBoolValue(v, config.temperature.bounded);
    }
    v = jsonValueAfterKey(block, "boundMode");
    if (v.empty()) v = jsonValueAfterKey(block, "boundsMode");
    if (!v.empty()) config.temperature.boundMode = valueUnquote(v);
    v = jsonValueAfterKey(block, "lowerBound");
    if (v.empty()) v = jsonValueAfterKey(block, "min");
    if (!v.empty()) {
        config.temperature.lowerBound =
            valueDoubleValue(v, config.temperature.lowerBound);
    }
    v = jsonValueAfterKey(block, "upperBound");
    if (v.empty()) v = jsonValueAfterKey(block, "max");
    if (!v.empty()) {
        config.temperature.upperBound =
            valueDoubleValue(v, config.temperature.upperBound);
    }
    v = jsonValueAfterKey(block, "bounds");
    const std::vector<double> bounds = valueNumbers(v);
    if (bounds.size() >= 2) {
        config.temperature.lowerBound = bounds[0];
        config.temperature.upperBound = bounds[1];
    }
};
parseTemperature(multiPhase);
parseTemperature(jsonBlockAfterKey(text, "temperature"));

auto parsePhaseChange = [&](const P* block) {
    if (!block) return;
    config.phaseChange.enabled = true;
    std::string v = jsonValueAfterKey(block, "enabled");
    if (v.empty()) v = jsonValueAfterKey(block, "enable");
    if (!v.empty()) {
        config.phaseChange.enabled =
            valueBoolValue(v, config.phaseChange.enabled);
    }
    v = jsonValueAfterKey(block, "model");
    if (v.empty()) v = jsonValueAfterKey(block, "type");
    if (!v.empty()) {
        config.phaseChange.model = valueUnquote(v);
    }
    v = jsonValueAfterKey(block, "saturationTemperature");
    if (v.empty()) v = jsonValueAfterKey(block, "Tsat");
    if (!v.empty()) {
        config.phaseChange.saturationTemperature =
            valueDoubleValue(v, config.phaseChange.saturationTemperature);
    }
    v = jsonValueAfterKey(block, "evaporationCoefficient");
    if (v.empty()) v = jsonValueAfterKey(block, "evaporationRate");
    if (v.empty()) v = jsonValueAfterKey(block, "coefficient");
    if (!v.empty()) {
        config.phaseChange.evaporationCoefficient =
            valueDoubleValue(v, config.phaseChange.evaporationCoefficient);
    }
    v = jsonValueAfterKey(block, "condensationCoefficient");
    if (v.empty()) v = jsonValueAfterKey(block, "condensationRate");
    if (v.empty()) v = jsonValueAfterKey(block, "coefficient");
    if (!v.empty()) {
        config.phaseChange.condensationCoefficient =
            valueDoubleValue(v, config.phaseChange.condensationCoefficient);
    }
    v = jsonValueAfterKey(block, "latentHeat");
    if (v.empty()) v = jsonValueAfterKey(block, "L");
    if (!v.empty()) {
        config.phaseChange.latentHeat =
            valueDoubleValue(v, config.phaseChange.latentHeat);
    }
    v = jsonValueAfterKey(block, "energyCoupling");
    if (!v.empty()) {
        config.phaseChange.energyCoupling =
            valueBoolValue(v, config.phaseChange.energyCoupling);
    }
};
parsePhaseChange(jsonBlockAfterKey(text, "phaseChange"));

}



void appendPhaseBlock(
        Config& config,
        const std::string& name,
        const P* block) {
    namespace MP = Physics::Multiphase;
    MP::PhaseProperties phase;
    phase.name = name;
    std::string v = jsonValueAfterKey(block, "model");
    if (v.empty()) v = jsonValueAfterKey(block, "thermoModel");
    if (v.empty()) v = jsonValueAfterKey(block, "equationOfState");
    if (v.empty()) v = jsonValueAfterKey(block, "eos");
    if (v.empty()) {
        const std::string typeText = jsonValueAfterKey(block, "type");
        const std::string typeKey = sourceTokenKey(typeText);
        if (typeKey == "perfectliquid" || typeKey == "perfectgas") {
            v = typeText;
        }
    }
    if (!v.empty()) {
        phase.thermoModel = valueUnquote(v);
        const std::string modelKey = sourceTokenKey(phase.thermoModel);
        if (modelKey == "perfectliquid") {
            phase.role = MP::PhaseRole::Liquid;
            phase.roleDeclared = true;
        } else if (modelKey == "perfectgas") {
            phase.role = MP::PhaseRole::Gas;
            phase.roleDeclared = true;
        }
    }
    v = jsonValueAfterKey(block, "role");
    if (v.empty()) {
        const std::string typeText = jsonValueAfterKey(block, "type");
        const std::string typeKey = sourceTokenKey(typeText);
        if (typeKey != "perfectliquid" && typeKey != "perfectgas") {
            v = typeText;
        }
    }
    if (!v.empty()) {
        phase.role = MP::parsePhaseRole(valueUnquote(v));
        phase.roleDeclared = phase.role != MP::PhaseRole::Unknown;
    }
    if (phase.role == MP::PhaseRole::Unknown) {
        phase.role = MP::parsePhaseRole(name);
    }
    v = jsonValueAfterKey(block, "density");
    if (v.empty()) v = jsonValueAfterKey(block, "rho0");
    if (v.empty()) v = jsonValueAfterKey(block, "rho");
    if (!v.empty()) phase.density = valueDoubleValue(v, phase.density);
    v = jsonValueAfterKey(block, "viscosity");
    if (v.empty()) v = jsonValueAfterKey(block, "mu");
    if (!v.empty()) {
        phase.viscosity = valueDoubleValue(v, phase.viscosity);
    }
    v = jsonValueAfterKey(block, "specificHeat");
    if (v.empty()) v = jsonValueAfterKey(block, "heatCapacity");
    if (v.empty()) v = jsonValueAfterKey(block, "Cp");
    if (v.empty()) v = jsonValueAfterKey(block, "cp");
    if (v.empty()) v = jsonValueAfterKey(block, "Cv");
    if (!v.empty()) {
        phase.specificHeat =
            valueDoubleValue(v, phase.specificHeat);
    }
    v = jsonValueAfterKey(block, "Cv");
    if (v.empty()) v = jsonValueAfterKey(block, "cv");
    if (!v.empty()) phase.Cv = valueDoubleValue(v, phase.Cv);
    v = jsonValueAfterKey(block, "e0");
    if (!v.empty()) phase.e0 = valueDoubleValue(v, phase.e0);
    v = jsonValueAfterKey(block, "T_ref");
    if (v.empty()) v = jsonValueAfterKey(block, "Tref");
    if (!v.empty()) phase.T_ref = valueDoubleValue(v, phase.T_ref);
    v = jsonValueAfterKey(block, "gamma");
    if (!v.empty()) phase.gamma = valueDoubleValue(v, phase.gamma);
    v = jsonValueAfterKey(block, "gasConstant");
    if (v.empty()) v = jsonValueAfterKey(block, "R");
    if (!v.empty()) phase.gasConstant = valueDoubleValue(v, phase.gasConstant);
    v = jsonValueAfterKey(block, "pInfinity");
    if (v.empty()) v = jsonValueAfterKey(block, "pInf");
    if (!v.empty()) phase.pInfinity = valueDoubleValue(v, phase.pInfinity);
    v = jsonValueAfterKey(block, "volumeFraction");
    if (v.empty()) v = jsonValueAfterKey(block, "alpha");
    if (!v.empty()) {
        phase.volumeFraction = valueDoubleValue(v, phase.volumeFraction);
        phase.initialAlphaDeclared = true;
    }
    v = jsonValueAfterKey(block, "temperature");
    if (v.empty()) v = jsonValueAfterKey(block, "T0");
    if (v.empty()) v = jsonValueAfterKey(block, "T");
    if (!v.empty()) {
        phase.initialTemperature = valueDoubleValue(v, phase.initialTemperature);
        phase.initialTemperatureDeclared = true;
    }
    v = jsonValueAfterKey(block, "velocity");
    if (v.empty()) v = jsonValueAfterKey(block, "U0");
    if (v.empty()) v = jsonValueAfterKey(block, "U");
    if (!v.empty()) {
        if (!valueVectorValue(v, &phase.initialVelocity)) {
            fatalCaseConfig("phase '" + name + "' has invalid velocity/U0.");
        }
        phase.initialVelocityDeclared = true;
    }
    phase.initialDensityDeclared = phase.density > 0.0;
    v = jsonValueAfterKey(block, "thermalConductivity");
    if (v.empty()) v = jsonValueAfterKey(block, "conductivity");
    if (v.empty()) v = jsonValueAfterKey(block, "lambda");
    if (v.empty()) v = jsonValueAfterKey(block, "k");
    if (!v.empty()) {
        phase.thermalConductivity =
            valueDoubleValue(v, phase.thermalConductivity);
    }
    auto parseTemperatureLaw = [&](const std::string& key,
                                   MP::TemperaturePropertyLaw& law) {
        const P* properties = jsonBlockAfterKey(block, "thermophysicalProperties");
        const P* lawBlock = jsonBlockAfterKey(properties, key);
        if (!lawBlock) lawBlock = jsonBlockAfterKey(block, key);
        if (!lawBlock) return;
        std::string lawValue = jsonValueAfterKey(lawBlock, "type");
        if (lawValue.empty()) {
            lawValue = jsonValueAfterKey(lawBlock, "model");
        }
        if (!lawValue.empty()) law.model = valueUnquote(lawValue);
        lawValue = jsonValueAfterKey(
            lawBlock, "referenceTemperature");
        if (lawValue.empty()) {
            lawValue = jsonValueAfterKey(lawBlock, "Tref");
        }
        if (!lawValue.empty()) {
            law.referenceTemperature = valueDoubleValue(
                lawValue, law.referenceTemperature);
        }
        lawValue = jsonValueAfterKey(lawBlock, "validTemperature");
        if (lawValue.empty()) {
            lawValue = jsonValueAfterKey(lawBlock, "temperatureRange");
        }
        const auto range = valueNumbers(lawValue);
        if (range.size() == 2) {
            law.minimumTemperature = range[0];
            law.maximumTemperature = range[1];
        }
        lawValue = jsonValueAfterKey(lawBlock, "coefficients");
        if (lawValue.empty()) {
            lawValue = jsonValueAfterKey(lawBlock, "coeffs");
        }
        law.coefficients = valueNumbers(lawValue);
    };
    parseTemperatureLaw("density", phase.densityLaw);
    parseTemperatureLaw("rho", phase.densityLaw);
    parseTemperatureLaw("dynamicViscosity", phase.viscosityLaw);
    parseTemperatureLaw("viscosity", phase.viscosityLaw);
    parseTemperatureLaw("mu", phase.viscosityLaw);
    parseTemperatureLaw("specificHeat", phase.specificHeatLaw);
    parseTemperatureLaw("Cp", phase.specificHeatLaw);
    parseTemperatureLaw(
        "thermalConductivity", phase.thermalConductivityLaw);
    parseTemperatureLaw("k", phase.thermalConductivityLaw);
    v = jsonValueAfterKey(block, "diameterModel");
    if (!v.empty()) phase.diameterModel = valueUnquote(v);
    const P* diameterBlock = jsonBlockAfterKey(block, phase.diameterModel + "Coeffs");
    if (!diameterBlock) {
        diameterBlock = jsonBlockAfterKey(block, "diameter");
    }
    v = jsonValueAfterKey(diameterBlock, "d");
    if (v.empty()) v = jsonValueAfterKey(diameterBlock, "diameter");
    if (!v.empty()) phase.diameter = valueDoubleValue(v, phase.diameter);
    v = jsonValueAfterKey(block, "residualAlpha");
    if (!v.empty()) {
        phase.residualAlpha =
            valueDoubleValue(v, phase.residualAlpha);
    }
    const P* speciesBlock = jsonBlockAfterKey(block, "species");
    for (const auto& [speciesName, speciesText] : jsonChildBlocks(speciesBlock)) {
        MP::SpeciesProperties species;
        species.name = speciesName;
        std::string y = jsonValueAfterKey(speciesText, "massFraction");
        if (y.empty()) y = jsonValueAfterKey(speciesText, "Y");
        if (!y.empty()) species.massFraction = valueDoubleValue(y, species.massFraction);
        phase.species.push_back(std::move(species));
    }
    config.phases.push_back(phase);
}

void parsePhaseDefinitions(
        const P* text,
        Config& config) {
    namespace MP = Physics::Multiphase;
const P* phases = jsonBlockAfterKey(text, "phases");
for (const auto& [name, block] : jsonChildBlocks(phases)) {
    appendPhaseBlock(config, name, block);
}

for (const std::string& phaseName : config.mixture.phaseNames) {
    const P* block = jsonBlockAfterKey(text, phaseName);
    if (block
        && MP::findPhase(config, phaseName) == nullptr) {
        appendPhaseBlock(config, phaseName, block);
    }
}


}



void parseEulerianProperties(
        const P* text,
        const P* multiPhase,
        Config& config) {
    namespace MP = Physics::Multiphase;
    std::string value;
if (MP::isEulerianEulerianType(config.type)) {
    value = jsonValueAfterKey(multiPhase, "phases");
    if (value.empty()) value = jsonValueAfterKey(text, "phases");
    config.eulerianEulerian.phaseNames = valueWordList(value);
    for (const std::string& phaseName :
         config.eulerianEulerian.phaseNames) {
        const P* block = jsonBlockAfterKey(text, phaseName);
        if (block && MP::findPhase(config, phaseName) == nullptr) {
            appendPhaseBlock(config, phaseName, block);
        }
    }
    value = jsonValueAfterKey(multiPhase, "referencePhase");
    if (value.empty()) value = jsonValueAfterKey(text, "referencePhase");
    config.eulerianEulerian.referencePhase = valueUnquote(value);
    value = jsonValueAfterKey(multiPhase, "geometryModel");
    if (value.empty()) value = jsonValueAfterKey(text, "geometryModel");
    const std::string geometryModel =
        FDM::normalizeToken(valueUnquote(value));
    if (!geometryModel.empty()
        && geometryModel != "cartesian"
        && geometryModel != "axisymmetric") {
        fatalCaseConfig(
            "eulerianEulerian geometryModel supports cartesian or "
            "axisymmetric.");
    }
    config.eulerianEulerian.axisymmetric =
        geometryModel == "axisymmetric";
    value = jsonValueAfterKey(multiPhase, "radialCoordinate");
    if (value.empty()) {
        value = jsonValueAfterKey(text, "radialCoordinate");
    }
    const std::string radial =
        FDM::normalizeToken(valueUnquote(value));
    if (!radial.empty()) {
        if (radial == "x") config.eulerianEulerian.radialCoordinate = 0;
        else if (radial == "y") config.eulerianEulerian.radialCoordinate = 1;
        else if (radial == "z") config.eulerianEulerian.radialCoordinate = 2;
        else fatalCaseConfig(
            "eulerianEulerian radialCoordinate supports x, y or z.");
    }
    if (!jsonValueAfterKey(multiPhase, "solver").empty()
        || !jsonValueAfterKey(text, "solver").empty()) {
        fatalCaseConfig(
            "eulerianEulerian phaseProperties must not select a solver; "
            "move type/algorithm to solvers/algorithm.yaml.");
    }

    const P* models = jsonBlockAfterKey(text, "interphaseModels");
    const P* pairs = jsonBlockAfterKey(models, "phasePairs");
    if (!pairs) pairs = jsonBlockAfterKey(models, "pairs");
    if (!pairs) {
        fatalCaseConfig(
            "eulerianEulerian interphaseModels must contain a phasePairs "
            "registry; the legacy single continuousPhase/dispersedPhase "
            "layout is no longer accepted.");
    }
    config.eulerianEulerian.phasePairs.clear();
    for (const auto& [pairName, pairBlock] : jsonChildBlocks(pairs)) {
        MP::PhasePairModelOptions pair;
        pair.name = pairName;
        value = jsonValueAfterKey(pairBlock, "continuousPhase");
        pair.continuousPhase = valueUnquote(value);
        value = jsonValueAfterKey(pairBlock, "dispersedPhase");
        pair.dispersedPhase = valueUnquote(value);

        const P* drag = jsonBlockAfterKey(pairBlock, "drag");
        const P* lift = jsonBlockAfterKey(pairBlock, "lift");
        const P* virtualMass = jsonBlockAfterKey(pairBlock, "virtualMass");
        const P* heat = jsonBlockAfterKey(pairBlock, "heatTransfer");
        const P* wall = jsonBlockAfterKey(pairBlock, "wallLubrication");
        const P* dispersion = jsonBlockAfterKey(pairBlock, "turbulentDispersion");
        const P* surfaceTension = jsonBlockAfterKey(pairBlock, "surfaceTension");

        value = jsonValueAfterKey(drag, "model");
        if (value.empty()) value = jsonValueAfterKey(drag, "type");
        if (!value.empty()) pair.dragModel = valueUnquote(value);
        value = jsonValueAfterKey(drag, "particleDiameter");
        if (value.empty()) value = jsonValueAfterKey(drag, "diameter");
        if (!value.empty()) {
            pair.particleDiameter =
                valueDoubleValue(value, pair.particleDiameter);
        }
        value = jsonValueAfterKey(lift, "model");
        if (value.empty()) value = jsonValueAfterKey(lift, "type");
        if (!value.empty()) pair.liftModel = valueUnquote(value);
        value = jsonValueAfterKey(lift, "coefficient");
        if (!value.empty()) {
            pair.liftCoefficient =
                valueDoubleValue(value, pair.liftCoefficient);
        }
        value = jsonValueAfterKey(virtualMass, "model");
        if (value.empty()) {
            value = jsonValueAfterKey(virtualMass, "type");
        }
        if (!value.empty()) {
            pair.virtualMassModel = valueUnquote(value);
        }
        value = jsonValueAfterKey(virtualMass, "coefficient");
        if (!value.empty()) {
            pair.virtualMassCoefficient =
                valueDoubleValue(value, pair.virtualMassCoefficient);
        }
        value = jsonValueAfterKey(heat, "model");
        if (value.empty()) value = jsonValueAfterKey(heat, "type");
        if (!value.empty()) {
            pair.heatTransferModel = valueUnquote(value);
        }
        value = jsonValueAfterKey(wall, "model");
        if (!value.empty()) {
            pair.wallLubricationModel = valueUnquote(value);
        }
        value = jsonValueAfterKey(wall, "wallPatch");
        if (value.empty()) value = jsonValueAfterKey(wall, "patch");
        if (!value.empty()) pair.wallPatch = valueUnquote(value);
        value = jsonValueAfterKey(wall, "C1");
        if (!value.empty()) {
            pair.wallLubricationC1 =
                valueDoubleValue(value, pair.wallLubricationC1);
        }
        value = jsonValueAfterKey(wall, "C2");
        if (!value.empty()) {
            pair.wallLubricationC2 =
                valueDoubleValue(value, pair.wallLubricationC2);
        }
        value = jsonValueAfterKey(dispersion, "model");
        if (!value.empty()) {
            pair.turbulentDispersionModel = valueUnquote(value);
        }
        value = jsonValueAfterKey(dispersion, "C0");
        if (value.empty()) value = jsonValueAfterKey(
            dispersion, "coefficient");
        if (!value.empty()) {
            pair.turbulentDispersionCoefficient =
                valueDoubleValue(
                    value, pair.turbulentDispersionCoefficient);
        }
        value = jsonValueAfterKey(
            dispersion, "turbulentKinematicViscosity");
        if (value.empty()) value = jsonValueAfterKey(dispersion, "nuT");
        if (!value.empty()) {
            pair.turbulentKinematicViscosity =
                valueDoubleValue(
                    value, pair.turbulentKinematicViscosity);
        }
        value = jsonValueAfterKey(dispersion, "sigmaAlpha");
        if (value.empty()) {
            value = jsonValueAfterKey(dispersion, "SchmidtNumber");
        }
        if (!value.empty()) {
            pair.turbulentSchmidtNumber =
                valueDoubleValue(value, pair.turbulentSchmidtNumber);
        }
        value = jsonValueAfterKey(surfaceTension, "model");
        if (value.empty()) {
            value = jsonValueAfterKey(surfaceTension, "type");
        }
        if (!value.empty()) {
            pair.surfaceTensionModel = valueUnquote(value);
        }
        value = jsonValueAfterKey(surfaceTension, "sigma");
        if (!value.empty()) {
            pair.surfaceTension =
                valueDoubleValue(value, pair.surfaceTension);
        }
        config.eulerianEulerian.phasePairs.push_back(std::move(pair));
    }

    const P* tensionRegistry = jsonBlockAfterKey(models, "surfaceTension");
    for (const auto& [name, tension] :
         jsonChildBlocks(tensionRegistry)) {
        const auto names =
            valueWordList(jsonValueAfterKey(tension, "phases"));
        if (names.size() != 2) {
            fatalCaseConfig(
                "surfaceTension entry '" + name
                + "' requires phases (phaseA phaseB).");
        }
        MP::PhasePairModelOptions* match = nullptr;
        for (auto& pair : config.eulerianEulerian.phasePairs) {
            const bool same =
                (MP::normalizePhaseName(pair.continuousPhase)
                     == MP::normalizePhaseName(names[0])
                 && MP::normalizePhaseName(pair.dispersedPhase)
                     == MP::normalizePhaseName(names[1]))
                || (MP::normalizePhaseName(pair.continuousPhase)
                        == MP::normalizePhaseName(names[1])
                    && MP::normalizePhaseName(pair.dispersedPhase)
                        == MP::normalizePhaseName(names[0]));
            if (same) match = &pair;
        }
        if (!match) {
            fatalCaseConfig(
                "surfaceTension entry '" + name
                + "' does not match a registered phase pair.");
        }
        value = jsonValueAfterKey(tension, "type");
        if (value.empty()) value = jsonValueAfterKey(tension, "model");
        match->surfaceTensionModel = valueUnquote(value);
        value = jsonValueAfterKey(tension, "sigma");
        match->surfaceTension =
            valueDoubleValue(value, match->surfaceTension);
    }
}


}



void parseSetAssignments(
        const P* text,
        Config& config) {
    namespace MP = Physics::Multiphase;
    std::string value;
for (const auto& [name, block] : jsonChildBlocks(text)) {
    const std::string prefix = "alpha.";
    if (name.rfind(prefix, 0) != 0 || name.size() <= prefix.size()) {
        continue;
    }
    appendPhaseBlock(config, name.substr(prefix.size()), block);
}

const P* sets = jsonBlockAfterKey(text, "sets");
for (const auto& [name, block] : jsonChildBlocks(sets)) {
    MP::SetPhaseAssignment assignment;
    assignment.setName = name;
    value = jsonValueAfterKey(block, "phase");
    if (value.empty()) value = jsonValueAfterKey(block, "value");
    assignment.phaseName = valueUnquote(value);
    config.setPhases.push_back(assignment);
}
if (sets != nullptr && sets->is_object()) {
    for (auto it = sets->begin(); it != sets->end(); ++it) {
        if (it.value().is_object()) continue;
        const std::string setName = it.key();
        std::size_t p = 0;
        const std::string tokenText = token(it.value());
        const std::string phaseName = readValueWord(tokenText, &p);
        if (!phaseName.empty()) {
            MP::SetPhaseAssignment assignment;
            assignment.setName = setName;
            assignment.phaseName = phaseName;
            config.setPhases.push_back(assignment);
        }
    }
}


}

} // namespace SF::IOPrivate::PhasePropertiesReader


namespace SF {

using namespace IOPrivate;


void CaseAdapter::decodeRuntimeSection() {
    const std::string context = "solvers/runtime.yaml";
    const P* text = &sections_.runtime;

    std::string value = jsonValueAfterKey(text, "createMesh");
    caseConfig_.createMesh = !value.empty() && valueBoolValue(value, false);

    value = jsonValueAfterKey(text, "startTime");
    if (!value.empty()) caseConfig_.time.startTime =
        valueDoubleValue(value, caseConfig_.time.startTime);
    value = jsonValueAfterKey(text, "endTime");
    if (!value.empty()) caseConfig_.time.endTime =
        valueDoubleValue(value, caseConfig_.time.endTime);
    value = jsonValueAfterKey(text, "endStep");
    if (!value.empty()) caseConfig_.time.endStep = std::max(
        0, valueIntValue(value, caseConfig_.time.endStep));
    // CFL / maxDeltaT 是数值 HOW，直接写进 typed numerics，不再经过
    // parser 全局暂存变量。
    value = jsonValueAfterKey(text, "CFL");
    if (!value.empty()) {
        caseConfig_.solver.numerics.cfl = valueDoubleValue(
            value, caseConfig_.solver.numerics.cfl);
    }
    value = jsonValueAfterKey(text, "maxDeltaT");
    if (!value.empty()) {
        caseConfig_.solver.numerics.maxDeltaT = valueDoubleValue(
            value, caseConfig_.solver.numerics.maxDeltaT);
    }

    value = jsonValueAfterKey(text, "writeInitial");
    if (!value.empty()) caseConfig_.writeInitial =
        valueBoolValue(value, caseConfig_.writeInitial);

    const std::string writeControl =
        valueLower(valueUnquote(jsonValueAfterKey(text, "writeControl")));
    if (writeControl != "timestep" && writeControl != "runtime") {
        fatalCaseConfig(
            context + " has unsupported writeControl '"
            + valueUnquote(jsonValueAfterKey(text, "writeControl"))
            + "'. Supported values are timeStep and runTime.");
    }
    value = jsonValueAfterKey(text, "writeInterval");
    if (!value.empty()) {
        if (writeControl == "timestep") {
            caseConfig_.time.writeByStep = true;
            caseConfig_.time.writeIntervalSteps = std::max(
                1, valueIntValue(
                    value, caseConfig_.time.writeIntervalSteps));
        } else {
            caseConfig_.time.writeByStep = false;
            caseConfig_.time.writeIntervalTime = valueDoubleValue(
                value, caseConfig_.time.writeIntervalTime);
        }
    }
}

void CaseAdapter::decodeOutputSection() {
    const std::string context = "solvers/runtime.yaml output";
    if (!sections_.hasOutput) {
        fatalCaseConfig(
            "solvers/runtime.yaml must declare output with jobName and "
            "outputDir.");
    }
    const P* text = &sections_.output;

    std::string value = jsonValueAfterKey(text, "jobName");
    if (value.empty()) {
        fatalCaseConfig(context + " must set jobName.");
    }
    jobName_ = valueUnquote(value);

    value = jsonValueAfterKey(text, "outputDir");
    if (value.empty()) {
        fatalCaseConfig(context + " must set outputDir.");
    }
    outputDir_ = valueUnquote(value);
}

void CaseAdapter::decodeNumericsSection() {
    const P* text = &sections_.numerics;
    FDM::NumericsConfig& numerics = caseConfig_.solver.numerics;

    std::string value;
    std::string convectionRecipeName;
    std::string diffusionRecipeName;

    // term recipe 是离散策略的名字，不是 scheme 字典的键。
    const P* terms = jsonBlockAfterKey(text, "terms");
    if (terms) {
        numerics.termRecipesDeclared = true;
        value = jsonValueAfterKey(terms, "convection");
        if (!value.empty()) convectionRecipeName = valueUnquote(value);
        value = jsonValueAfterKey(terms, "diffusion");
        if (!value.empty()) diffusionRecipeName = valueUnquote(value);
    }
    if (numerics.termRecipesDeclared) {
        // Native term 输入是完整选择边界：缺失的 role 必须保持缺失，让
        // NumericalCompiler 拒绝没有 recipe 的数学 term，而不是继承隐藏默认。
        numerics.recipes.convection.reset();
        numerics.recipes.diffusion.reset();
        if (!FDM::normalizeToken(convectionRecipeName).empty()) {
            numerics.recipes.convection =
                FDM::resolveConvectionTermRecipe(convectionRecipeName);
            numerics.convection = numerics.recipes.convection->convection();
            numerics.reconstruction =
                numerics.recipes.convection->reconstruction();
            numerics.flux = numerics.recipes.convection->flux();
            if (numerics.recipes.convection->id()
                == FDM::TermRecipeId::PrimitiveUpwind1) {
                numerics.formulation =
                    FDM::EquationFormulation::PrimitiveDifferential;
            }
        }
        if (!FDM::normalizeToken(diffusionRecipeName).empty()) {
            numerics.recipes.diffusion =
                FDM::resolveDiffusionTermRecipe(diffusionRecipeName);
            numerics.viscous = numerics.recipes.diffusion->diffusion();
        }
    }

    const P* time = jsonBlockAfterKey(text, "time");
    if (time) {
        value = jsonValueAfterKey(time, "default");
        if (!value.empty()) {
            numerics.timeRecipe =
                FDM::resolveTimeRecipe(valueUnquote(value));
            numerics.timeRecipeDeclared = true;
        }
        numerics.recipes.time = numerics.timeRecipe;
    }

    const P* convection = jsonBlockAfterKey(text, "convection");
    if (convection) {
        value = jsonValueAfterKey(convection, "default");
        if (!value.empty()) {
            numerics.convection =
                FDM::parseConvectionScheme(valueUnquote(value));
        }
        value = jsonValueAfterKey(convection, "formulation");
        if (!value.empty()) {
            numerics.formulation =
                FDM::parseEquationFormulation(valueUnquote(value));
        }
        value = jsonValueAfterKey(convection, "reconstruction");
        if (!value.empty()) {
            numerics.reconstruction =
                FDM::parseReconstructionVariable(valueUnquote(value));
        }
        value = jsonValueAfterKey(convection, "flux");
        if (!value.empty()) {
            numerics.flux = FDM::parseFluxSplitter(valueUnquote(value));
        }
        value = jsonValueAfterKey(convection, "interfaceFlux");
        if (!value.empty()) {
            numerics.interfaceFlux =
                FDM::parseInterfaceFluxPolicy(valueUnquote(value));
        }
    }

    const P* diffusion = jsonBlockAfterKey(text, "diffusion");
    if (diffusion) {
        value = jsonValueAfterKey(diffusion, "default");
        if (!value.empty()) {
            numerics.viscous =
                FDM::parseViscousScheme(valueUnquote(value));
        }
    }

    const P* sources = jsonBlockAfterKey(text, "sources");
    if (sources) {
        value = jsonValueAfterKey(sources, "default");
        if (!value.empty()) sourceScheme_ = valueUnquote(value);
    }

    // transport / pressureCorrection 是数值控制语义块。旧的外部字典
    // viscous / pressureBased 只是它们的字典形状，不再是解码入口。
    const P* transport = jsonBlockAfterKey(text, "transport");
    if (transport) {
        value = jsonValueAfterKey(transport, "type");
        if (!value.empty()) {
            numerics.viscous =
                FDM::parseViscousScheme(valueUnquote(value));
        }
        value = jsonValueAfterKey(transport, "enable");
        if (!value.empty()) {
            numerics.viscousEnabled =
                valueBoolValue(value, numerics.viscousEnabled);
        }
        value = jsonValueAfterKey(transport, "mu");
        if (!value.empty()) {
            numerics.dynamicViscosity = valueDoubleValue(
                value, numerics.dynamicViscosity);
        }
        value = jsonValueAfterKey(transport, "Pr");
        if (!value.empty()) {
            numerics.prandtl = valueDoubleValue(value, numerics.prandtl);
        }
    }

    const P* pressure = jsonBlockAfterKey(text, "pressureCorrection");
    if (pressure) {
        FDM::PressureCorrectionConfig& pressureConfig =
            caseConfig_.solver.pressure;
        value = jsonValueAfterKey(pressure, "maxIterations");
        if (!value.empty()) pressureConfig.maxIterations =
            valueIntValue(value, pressureConfig.maxIterations);
        value = jsonValueAfterKey(pressure, "relativeTolerance");
        if (!value.empty()) pressureConfig.relativeTolerance =
            valueDoubleValue(value, pressureConfig.relativeTolerance);
        value = jsonValueAfterKey(pressure, "absoluteTolerance");
        if (!value.empty()) pressureConfig.absoluteTolerance =
            valueDoubleValue(value, pressureConfig.absoluteTolerance);
        value = jsonValueAfterKey(pressure, "pressureRelaxation");
        if (!value.empty()) pressureConfig.relaxation =
            valueDoubleValue(value, pressureConfig.relaxation);
        value = jsonValueAfterKey(pressure, "velocityRelaxation");
        if (!value.empty()) pressureConfig.velocityRelaxation =
            valueDoubleValue(value, pressureConfig.velocityRelaxation);
    }

    // source 启用状态只能由 sourceScheme 推导；数值块不选择物理模型。
    // 最终选择由 buildCaseConfig() 转成 typed SourceKind 列表。
}

void CaseAdapter::decodeParallelSection() {
    if (!sections_.hasParallel) return;
    const P* parallel = &sections_.parallel;

    std::string value = jsonValueAfterKey(parallel, "enabled");
    if (!value.empty()) caseConfig_.parallel.enabled = valueBoolValue(
        value, caseConfig_.parallel.enabled);
    value = jsonValueAfterKey(parallel, "nProcs");
    if (!value.empty()) caseConfig_.parallel.processCount = valueIntValue(
        value, caseConfig_.parallel.processCount);
    value = jsonValueAfterKey(parallel, "split");
    std::array<int, 3> split{};
    if (!value.empty() && valueIntegerTriple(value, &split)) {
        caseConfig_.parallel.partitionSplit = split;
        caseConfig_.parallel.automaticPartition = false;
        caseConfig_.parallel.processCount = split[0] * split[1] * split[2];
    } else if (!value.empty() && valueLower(valueUnquote(value)) == "auto") {
        caseConfig_.parallel.automaticPartition = true;
    }
    value = jsonValueAfterKey(parallel, "haloTolerance");
    if (!value.empty()) {
        caseConfig_.parallel.haloTolerance = valueDoubleValue(
            value, caseConfig_.parallel.haloTolerance);
    }
}

void CaseAdapter::decodeAlgorithmSection() {
    const std::string context = "solvers/algorithm.yaml";
    if (!sections_.hasAlgorithm) {
        fatalCaseConfig(
            "solvers/solvers.yaml must declare an algorithm entry backed by "
            "solvers/algorithm.yaml with type densityBase or pressureBase.");
    }
    const P* text = &sections_.algorithm;
    FDM::PressureCorrectionConfig config;
    std::string value = jsonValueAfterKey(text, "type");
    const std::string type = FDM::normalizeToken(valueUnquote(value));
    // Legacy 兼容标签：只在这里翻译，向 composition root 传递"请求哪一组方程"。
    // 它不再是 runtime 求解器身份，也不进入 ResolvedSimulationSystem。
    if (!sections_.nativeAlgorithm.empty()) {
        legacyFlowLabel_.clear();
    } else if (type == "pressurebase") {
        legacyFlowLabel_ = "pressureBase";
    } else if (type == "densitybase") {
        legacyFlowLabel_ = "densityBase";
    } else {
        fatalCaseConfig(
            context + " requires type densityBase or pressureBase.");
    }

    value = jsonValueAfterKey(text, "algorithm");
    const std::string algorithm = FDM::normalizeToken(valueUnquote(value));
    if (algorithm == "simple") config.coupling.preset = FDM::PressureCouplingPreset::SIMPLE;
    else if (algorithm == "piso") config.coupling.preset = FDM::PressureCouplingPreset::PISO;
    else if (algorithm == "pimple") config.coupling.preset = FDM::PressureCouplingPreset::PIMPLE;
    else if ((!algorithm.empty() && algorithm!="explicit") || legacyFlowLabel_ == "pressureBase") {
        fatalCaseConfig(
            context + " requires a known coupling algorithm "
            "SIMPLE, PISO, or PIMPLE.");
    }

    caseConfig_.pressureCouplingDeclared=!algorithm.empty() && algorithm!="explicit";
    const P* algorithmBlock = jsonBlockAfterKey(text, algorithm.empty() ? "PIMPLE" : valueUnquote(value));
    const P* controls = algorithmBlock ? algorithmBlock : text;
    auto readInt = [&](const char* key, int& target) {
        const std::string raw = jsonValueAfterKey(controls, key);
        if (!raw.empty()) target = valueIntValue(raw, target);
    };
    auto readDouble = [&](const char* key, double& target) {
        const std::string raw = jsonValueAfterKey(controls, key);
        if (!raw.empty()) target = valueDoubleValue(raw, target);
    };
    readInt("outerCorrectors", config.coupling.outerCorrectors);
    readInt("pressureCorrectors", config.coupling.pressureCorrectors);
    readInt("nonOrthogonalCorrectors", config.coupling.nonOrthogonalCorrectors);
    readDouble("momentumRelaxation", config.coupling.momentumRelaxation);
    readDouble("pressureRelaxation", config.coupling.pressureRelaxation);
    readDouble("relativeTolerance", config.relativeTolerance);
    readDouble("absoluteTolerance", config.absoluteTolerance);
    readDouble("phaseSourceCFL", config.phaseTransport.sourceCfl);
    const std::string phaseConvection = FDM::normalizeToken(
        valueUnquote(jsonValueAfterKey(controls, "phaseConvection")));
    if (!phaseConvection.empty() && phaseConvection != "upwind") {
        fatalCaseConfig(
            context + " phaseConvection currently supports only upwind.");
    }
    readInt("referenceCell", config.reference.referenceCell);
    readDouble("referencePressure", config.reference.referencePressure);

    const P* linear = jsonBlockAfterKey(text, "linearSolvers");
    auto parseLinear = [&](const char* name, FDM::LinearSolverConfig& target) {
        const P* block = jsonBlockAfterKey(linear, name);
        if (!block) {
            fatalCaseConfig(
                context + " linearSolvers is missing '"
                + name + "'.");
        }
        std::string raw = jsonValueAfterKey(block, "backend");
        if (FDM::normalizeToken(valueUnquote(raw)) != "hypre") {
            fatalCaseConfig(
                std::string("linear solver '") + name
                + "' requires backend hypre.");
        }
        raw = jsonValueAfterKey(block, "solver");
        const std::string method = FDM::normalizeToken(valueUnquote(raw));
        if (method == "flexgmres" || method == "fgmres") {
            target.method = FDM::KrylovMethod::FlexGMRES;
        } else if (method == "pcg" || method == "cg") {
            target.method = FDM::KrylovMethod::PCG;
        } else {
            fatalCaseConfig(
                std::string("linear solver '") + name
                + "' supports flexGMRES or PCG.");
        }
        raw = jsonValueAfterKey(block, "preconditioner");
        const std::string preconditioner =
            FDM::normalizeToken(valueUnquote(raw));
        if (preconditioner == "boomeramg") {
            target.preconditioner=FDM::LinearPreconditioner::BoomerAMG;
        } else if (preconditioner == "ilu") {
            target.preconditioner=FDM::LinearPreconditioner::ILU;
        } else if (preconditioner == "kktblockschur") {
            target.preconditioner=FDM::LinearPreconditioner::KKTBlockSchur;
        } else {
            fatalCaseConfig(
                std::string("linear solver '") + name
                + "' preconditioner supports boomerAMG, ILU or kktBlockSchur.");
        }
        raw = jsonValueAfterKey(block, "equilibration");
        const std::string equilibration =
            FDM::normalizeToken(valueUnquote(raw));
        if (equilibration.empty() || equilibration == "none") {
            target.equilibration = FDM::LinearEquilibration::None;
        } else if (equilibration == "rowmax") {
            target.equilibration = FDM::LinearEquilibration::RowMax;
        } else {
            fatalCaseConfig(
                std::string("linear solver '") + name
                + "' equilibration supports none or rowMax.");
        }
        auto blockInt = [&](const char* key, int& destination) {
            const std::string item = jsonValueAfterKey(block, key);
            if (!item.empty()) destination = valueIntValue(item, destination);
        };
        auto blockDouble = [&](const char* key, double& destination) {
            const std::string item = jsonValueAfterKey(block, key);
            if (!item.empty()) destination = valueDoubleValue(item, destination);
        };
        blockInt("maxIterations", target.maxIterations);
        blockInt("krylovDimension", target.krylovDimension);
        blockInt("structureRebuildInterval",
                 target.structureRebuildInterval);
        blockInt("rebuildInterval", target.structureRebuildInterval);
        blockInt("preconditionerRefreshInterval",
                 target.preconditionerRefreshInterval);
        if(target.preconditioner==FDM::LinearPreconditioner::ILU) {
            const std::string type=jsonValueAfterKey(block,"iluType");
            const std::string fill=jsonValueAfterKey(block,"iluLevelOfFill");
            if(type.empty()||fill.empty())fatalCaseConfig(
                std::string("linear solver '")+name
                +"' with ILU requires explicit iluType and iluLevelOfFill.");
            target.iluType=valueIntValue(type,target.iluType);
            target.iluLevelOfFill=valueIntValue(fill,target.iluLevelOfFill);
        }
        blockDouble("relativeTolerance", target.relativeTolerance);
        blockDouble("absoluteTolerance", target.absoluteTolerance);
        blockInt("amgCoarsenType",target.amgCoarsenType);
        blockInt("amgRelaxType",target.amgRelaxType);
        blockInt("amgSweeps",target.amgSweeps);
        blockInt("amgMaxLevels",target.amgMaxLevels);
        blockDouble("amgStrongThreshold",target.amgStrongThreshold);
        blockInt("schurDenseLimit",target.schurDenseLimit);
    };
    if (legacyFlowLabel_=="pressureBase"
        || (!sections_.nativeAlgorithm.empty() && caseConfig_.pressureCouplingDeclared)) {
        parseLinear("pressure", config.linear.pressure);
        parseLinear("momentum", config.linear.momentum);
        parseLinear("energy", config.linear.energy);
        config.linear.turbulence = config.linear.energy;
        if (jsonBlockAfterKey(linear, "turbulence")) {
            parseLinear("turbulence", config.linear.turbulence);
        }
        try {
            FDM::validatePressureCorrectionConfig(config);
        } catch (const std::exception& e) {
            fatalCaseConfig(e.what());
        }
    }
    solverProperties_ = config;
    solverPropertiesLoaded_ = true;
}


void CaseAdapter::decodeFields() {
    // field 的存在性就是它的语义声明；不再从 <startTime>/ 目录名或文件
    // 后缀推断，也不存在第二条 field authority。
    if (const Model::FieldDescriptor* field = sections_.field("U")) {
        decodeVectorField(*field, caseConfig_.solver.initial.velocity,
                          caseConfig_.solver.boundaries.velocity, nullptr);
    }
    if (const Model::FieldDescriptor* field = sections_.field("p")) {
        decodeScalarField(*field, caseConfig_.solver.initial.pressure,
                          caseConfig_.solver.boundaries.pressure,
                          nullptr);
    }
    if (const Model::FieldDescriptor* field = sections_.field("rho")) {
        decodeScalarField(*field, caseConfig_.solver.initial.density,
                          caseConfig_.solver.boundaries.density, nullptr);
    }
    if (const Model::FieldDescriptor* field = sections_.field("rhoE")) {
        std::vector<BCSetting<double>> unusedBoundary;
        decodeScalarField(*field, caseConfig_.solver.initial.energy,
                          unusedBoundary, nullptr);
    }
    if (const Model::FieldDescriptor* field = sections_.field("T")) {
        decodeTemperatureField(*field);
    }
    if (const Model::FieldDescriptor* field = sections_.field("k")) {
        decodeScalarField(*field,
                          caseConfig_.solver.turbulence.scalars.kInitial,
                          caseConfig_.solver.turbulence.scalars.kBoundary,
                          nullptr);
    }
    if (const Model::FieldDescriptor* field = sections_.field("epsilon")) {
        decodeScalarField(
            *field,
            caseConfig_.solver.turbulence.scalars.epsilonInitial,
            caseConfig_.solver.turbulence.scalars.epsilonBoundary,
            nullptr);
    }
    if (const Model::FieldDescriptor* field = sections_.field("omega")) {
        decodeScalarField(
            *field,
            caseConfig_.solver.turbulence.scalars.omegaInitial,
            caseConfig_.solver.turbulence.scalars.omegaBoundary,
            nullptr);
    }
    if (const Model::FieldDescriptor* field = sections_.field("phi")) {
        decodePhiField(*field);
    }
    if (multiPhaseLoaded_ && Physics::Multiphase::isThermalType(multiPhaseConfig_.type)
        && multiPhaseConfig_.temperature.fieldName!="T") {
        if (const auto* field=sections_.field(multiPhaseConfig_.temperature.fieldName))
            decodeScalarField(*field,multiPhaseConfig_.temperature.initialConditions,
                multiPhaseConfig_.temperature.boundaryConditions,&multiPhaseConfig_.temperature.defaultValue);
    }
    bool parsedAlpha = false;
    if (const Model::FieldDescriptor* field = sections_.field("alpha")) {
        if (!multiPhaseLoaded_) {
            std::cerr << "[SF FATAL] fields/alpha requires explicit "
                      << "models/multiPhase.yaml with a phaseSystem. "
                      << "Do not enable a mixture model implicitly."
                      << std::endl;
            std::exit(1);
        }
        decodeScalarField(
            *field,
            multiPhaseConfig_.alpha.initialConditions,
            multiPhaseConfig_.alpha.boundaryConditions,
            &multiPhaseConfig_.alpha.defaultValue);
        parsedAlpha = true;
    }
    if (multiPhaseLoaded_
        && Physics::Multiphase::isMixtureType(multiPhaseConfig_.type)) {
        for (const std::string& fieldName :
             multiPhaseConfig_.mixture.phaseFractionFields) {
            if (parsedAlpha) break;
            const Model::FieldDescriptor* field = sections_.field(fieldName);
            if (!field) continue;
            multiPhaseConfig_.alpha.fieldName = fieldName;
            const size_t dot = fieldName.find_last_of('.');
            if (dot != std::string::npos && dot + 1 < fieldName.size()) {
                multiPhaseConfig_.alpha.phaseName = fieldName.substr(dot + 1);
            }
            decodeScalarField(
                *field,
                multiPhaseConfig_.alpha.initialConditions,
                multiPhaseConfig_.alpha.boundaryConditions,
                &multiPhaseConfig_.alpha.defaultValue);
            parsedAlpha = true;
        }
        if (!parsedAlpha) {
            fatalCaseConfig(
                "a mixture phaseSystem requires at least one phase-fraction "
                "field in fields/, e.g. alpha or alpha.liquid.");
        }
    }
    if (multiPhaseLoaded_
        && Physics::Multiphase::isEulerianEulerianType(
            multiPhaseConfig_.type)) {
        for (const std::string& phaseName :
             multiPhaseConfig_.eulerianEulerian.phaseNames) {
            auto phaseIt = std::find_if(
                multiPhaseConfig_.phases.begin(),
                multiPhaseConfig_.phases.end(),
                [&](const Physics::Multiphase::PhaseProperties& phase) {
                    return Physics::Multiphase::normalizePhaseName(phase.name)
                        == Physics::Multiphase::normalizePhaseName(phaseName);
                });
            if (phaseIt == multiPhaseConfig_.phases.end()) {
                fatalCaseConfig(
                    "Eulerian-Eulerian field references undeclared phase '"
                    + phaseName + "'.");
            }
            auto& phase = *phaseIt;
            const std::string suffix = "." + phaseName;
            const Model::FieldDescriptor* alpha = sections_.field("alpha" + suffix);
            const Model::FieldDescriptor* density = sections_.field("rho" + suffix);
            const Model::FieldDescriptor* velocity = sections_.field("U" + suffix);
            const Model::FieldDescriptor* temperature = sections_.field("T" + suffix);
            if (!alpha || !density || !velocity || !temperature) {
                fatalCaseConfig(
                    "Eulerian-Eulerian requires alpha.<phase>, rho.<phase>, "
                    "U.<phase>, and T.<phase> for phase '" + phaseName + "'.");
            }
            decodeScalarField(
                *alpha, phase.alphaInitial, phase.alphaBoundary,
                &phase.volumeFraction);
            decodeScalarField(
                *density, phase.densityInitial, phase.densityBoundary,
                &phase.density);
            decodeVectorField(
                *velocity, phase.velocityInitial, phase.velocityBoundary,
                &phase.initialVelocity);
            decodeScalarField(
                *temperature, phase.temperatureInitial,
                phase.temperatureBoundary, &phase.initialTemperature);
            phase.initialAlphaDeclared = true;
            phase.initialDensityDeclared = true;
            phase.initialVelocityDeclared = true;
            phase.initialTemperatureDeclared = true;
        }
        try {
            Physics::Multiphase::validateMultiPhaseConfig(
                multiPhaseConfig_, "Eulerian-Eulerian fields/");
        } catch (const std::exception& e) {
            fatalCaseConfig(e.what());
        }
    }
}


void CaseAdapter::decodeVectorField(
    const Model::FieldDescriptor& field,
    std::vector<BCSetting<Vector3>>& internal,
    std::vector<BCSetting<Vector3>>& boundary,
    Vector3* uniformValue) {
    std::vector<BCSetting<Vector3>> ic;
    std::vector<BCSetting<Vector3>> bc;

    Vector3 value;
    std::string raw = token(field.initial);
    if (!raw.empty() && valueVectorValue(raw, &value)) {
        if (uniformValue != nullptr) *uniformValue = value;
        ic.push_back({"all", FIXED_VALUE, value});
    }

    const P* sets = &field.sets;
    for (const auto& [name, block] : jsonChildBlocks(sets)) {
        raw = jsonValueAfterKey(block, "value");
        Vector3 setValue;
        if (!valueVectorValue(raw, &setValue)) setValue = Vector3(0, 0, 0);
        BCType type = valueBCTypeValue(jsonValueAfterKey(block, "type"),
                                      FIXED_VALUE);
        ic.push_back({name, type, setValue});
    }

    const P* patches = &field.boundaries;
    for (const auto& [name, block] : jsonChildBlocks(patches)) {
        raw = jsonValueAfterKey(block, "value");
        Vector3 patchValue;
        if (!valueVectorValue(raw, &patchValue)) patchValue = Vector3(0, 0, 0);
        const std::string wallType = valueLower(valueUnquote(jsonValueAfterKey(block,"type")));
        if (wallType == "noslip" || wallType == "slip") {
            if (!raw.empty() && (!valueVectorValue(raw,&patchValue)
                || patchValue.x!=0. || patchValue.y!=0. || patchValue.z!=0.))
                throw std::runtime_error(field.name+": boundary '"+name+"' uses stationary "+wallType
                    +" and cannot prescribe a nonzero or malformed wall velocity.");
            patchValue = {};
        }
        BCType type = valueBCTypeValue(jsonValueAfterKey(block, "type"),
                                      ZERO_GRADIENT,true);
        bc.push_back({name, type, patchValue});
    }

    if (!ic.empty()) internal = ic;
    if (!bc.empty()) boundary = bc;
}

void CaseAdapter::decodeScalarField(
    const Model::FieldDescriptor& field,
    std::vector<BCSetting<double>>& internal,
    std::vector<BCSetting<double>>& boundary,
    double* uniformValue) {
    std::vector<BCSetting<double>> ic;
    std::vector<BCSetting<double>> bc;

    double value = 0.0;
    std::string raw = token(field.initial);
    if (!raw.empty() && valueScalarValue(raw, &value)) {
        if (uniformValue != nullptr) *uniformValue = value;
        ic.push_back({"all", FIXED_VALUE, value});
    }

    const P* sets = &field.sets;
    for (const auto& [name, block] : jsonChildBlocks(sets)) {
        raw = jsonValueAfterKey(block, "value");
        double setValue = 0.0;
        valueScalarValue(raw, &setValue);
        BCType type = valueBCTypeValue(jsonValueAfterKey(block, "type"),
                                      FIXED_VALUE);
        ic.push_back({name, type, setValue});
    }

    const P* patches = &field.boundaries;
    for (const auto& [name, block] : jsonChildBlocks(patches)) {
        raw = jsonValueAfterKey(block, "value");
        double patchValue = 0.0;
        valueScalarValue(raw, &patchValue);
        BCType type = valueBCTypeValue(jsonValueAfterKey(block, "type"),
                                      ZERO_GRADIENT);
        bc.push_back({name, type, patchValue});
    }

    if (!ic.empty()) internal = ic;
    if (!bc.empty()) boundary = bc;
}

void CaseAdapter::decodeTemperatureField(const Model::FieldDescriptor& field) {
    std::vector<BCSetting<double>> ic;
    std::vector<ThermalBCSetting> bc;
    wallHeatBoundarySettings_.clear();

    double value = 0.0;
    std::string raw = token(field.initial);
    if (!raw.empty() && valueScalarValue(raw, &value)) {
        ic.push_back({"all", FIXED_VALUE, value});
        if (multiPhaseLoaded_) {
            multiPhaseConfig_.temperature.defaultValue = value;
        }
    }

    const P* sets = &field.sets;
    for (const auto& [name, block] : jsonChildBlocks(sets)) {
        raw = jsonValueAfterKey(block, "value");
        double setValue = 0.0;
        valueScalarValue(raw, &setValue);
        ic.push_back({name, FIXED_VALUE, setValue});
    }

    const P* patches = &field.boundaries;
    for (const auto& [name, block] : jsonChildBlocks(patches)) {
        const std::string typeText = jsonValueAfterKey(block, "type");
        if (valueWallHeatSourceType(typeText)) {
            WallHeatSetting setting;
            setting.patch = name;
            setting.mode = valueUnquote(typeText);
            raw = jsonValueAfterKey(block, "q");
            if (raw.empty()) raw = jsonValueAfterKey(block, "wallHeatFlux");
            if (raw.empty()) raw = jsonValueAfterKey(block, "heatFlux");
            if (raw.empty()) raw = jsonValueAfterKey(block, "value");
            if (raw.empty()) {
                fatalCaseConfig(field.name + ": patch '" + name
                                + "' uses wallHeatFlux but does not provide q.");
            }
            setting.heatFlux = valueDoubleValue(raw, setting.heatFlux);
            raw = jsonValueAfterKey(block, "coupleEnergy");
            if (raw.empty()) raw = jsonValueAfterKey(block, "energyCoupling");
            if (!raw.empty()) setting.coupleEnergy = valueBoolValue(raw, true);
            wallHeatBoundarySettings_.push_back(setting);
            // wallHeatFlux 是外部能量输入。RPI 只用它划分蒸发质量通量，
            // 不得因为相变模型消费了 q 就从总能量方程删除该边界功率。
            appendSourceSchemeToken(sourceScheme_, "WallHeat");
            bc.push_back({name, ThermalBCType::ZeroGradient, 0.0});
            continue;
        }

        const ThermalBCType type = valueThermalBCTypeValue(
            typeText,
            ThermalBCType::ZeroGradient);
        double patchValue = 0.0;
        if (type == ThermalBCType::HeatFlux) {
            raw = jsonValueAfterKey(block, "q");
            if (raw.empty()) raw = jsonValueAfterKey(block, "heatFlux");
            if (raw.empty()) raw = jsonValueAfterKey(block, "value");
        } else {
            raw = jsonValueAfterKey(block, "value");
        }
        valueScalarValue(raw, &patchValue);
        bc.push_back({name, type, patchValue});
    }

    std::vector<BCSetting<double>> scalarBC;
    scalarBC.reserve(bc.size());
    for (const ThermalBCSetting& condition : bc) {
        switch (condition.type) {
        case ThermalBCType::FixedTemperature:
            scalarBC.push_back({condition.name, FIXED_VALUE,
                                condition.value});
            break;
        case ThermalBCType::ZeroGradient:
        case ThermalBCType::Adiabatic:
            scalarBC.push_back({condition.name, ZERO_GRADIENT, 0.0});
            break;
        case ThermalBCType::Empty:
            scalarBC.push_back({condition.name, EMPTY, 0.0});
            break;
        case ThermalBCType::HeatFlux:
            break;
        }
    }
    temperatureOutputBC_ = scalarBC;
    if (multiPhaseLoaded_ && (!Physics::Multiphase::isThermalType(multiPhaseConfig_.type)
        || multiPhaseConfig_.temperature.fieldName=="T")) {
        multiPhaseConfig_.temperature.initialConditions = ic;
        multiPhaseConfig_.temperature.boundaryConditions = scalarBC;
        return;
    }

    if (!ic.empty()) caseConfig_.solver.initial.temperature = ic;
    if (!bc.empty()) caseConfig_.solver.boundaries.thermal = bc;
}


void CaseAdapter::decodeTurbulenceSection() {
    if (!sections_.hasTurbulence) return;
    const std::string context = "models/turbulence.yaml";
    const P* text = &sections_.turbulence;
    FDM::TurbulenceConfig& turbulence = caseConfig_.solver.turbulence;

    std::string value;
    turbulence.enabled = modelEnabled(text);
    if (!turbulence.enabled) {
        turbulence.family = FDM::TurbulenceFamily::DNS;
        turbulence.model = FDM::TurbulenceModelKind::DNS;
        return;
    }

    value = jsonValueAfterKey(text, "simulationType");
    if (value.empty()) value = jsonValueAfterKey(text, "type");
    if (value.empty()) {
        fatalCaseConfig(
            "models/turbulence.yaml has enabled=true but no "
            "simulationType/type. Use RAS, LES, or DNS.");
    }
    std::string family;
    std::string familyText;
    if (!value.empty()) {
        familyText = valueUnquote(value);
        family = upperToken(familyText);
        if (family == "LAMINAR" || family == "DNS") {
            turbulence.enabled = false;
            turbulence.family = FDM::TurbulenceFamily::DNS;
            turbulence.model = FDM::TurbulenceModelKind::DNS;
            return;
        }
        if (family != "RAS" && family != "LES") {
            fatalCaseConfig(
                "models/turbulence.yaml simulationType must be "
                "RAS, LES, or DNS; received '" + familyText + "'.");
        }
        turbulence.family = FDM::parseTurbulenceFamily(familyText);
    }

    const P* modelBlock = jsonBlockAfterKey(text,family);
    if (!modelBlock) {
        fatalCaseConfig(
            "models/turbulence.yaml selects "+familyText
            +" but has no matching "+family+" model block.");
    }
    if (modelBlock) {
        value = jsonValueAfterKey(modelBlock, "model");
        if (value.empty()) value = jsonValueAfterKey(modelBlock, "type");
        if (value.empty()) {
            fatalCaseConfig(
                "models/turbulence.yaml RAS/LES block must provide "
                "model/type.");
        }
        turbulence.model =
            FDM::parseTurbulenceModelKind(valueUnquote(value));

        auto setCoeff = [&](const char* key, double& target) {
            std::string coeff = jsonValueAfterKey(modelBlock, key);
            if (!coeff.empty()) target = valueDoubleValue(coeff, target);
        };
        FDM::TurbulenceCoefficients& coefficients = turbulence.coefficients;
        setCoeff("Cmu", coefficients.cMu);
        setCoeff("cMu", coefficients.cMu);
        setCoeff("C1", coefficients.c1);
        setCoeff("c1", coefficients.c1);
        setCoeff("C2", coefficients.c2);
        setCoeff("c2", coefficients.c2);
        setCoeff("sigmaK", coefficients.sigmaK);
        setCoeff("sigmaEpsilon", coefficients.sigmaEpsilon);
        setCoeff("betaStar", coefficients.betaStar);
        setCoeff("beta1", coefficients.beta1);
        setCoeff("gamma1", coefficients.gamma1);
        setCoeff("a1", coefficients.a1);
        setCoeff("sigmaOmega1", coefficients.sigmaOmega1);
        setCoeff("cSmagorinsky", coefficients.cSmagorinsky);
        setCoeff("filterScale", coefficients.filterScale);
        setCoeff("kFloor", coefficients.kFloor);
        setCoeff("epsilonFloor", coefficients.epsilonFloor);
        setCoeff("omegaFloor", coefficients.omegaFloor);
    }
    turbulence.phaseNames =
        valueWordList(jsonValueAfterKey(text, "phases"));
    turbulence.wallPatches =
        valueWordList(jsonValueAfterKey(text, "walls"));
    value=jsonValueAfterKey(text,"transportMethod");
    if(!value.empty()) {
        const auto method=valueUnquote(value);
        if(method!="stageConservative" && method!="legacyOnceStep")throw std::runtime_error("Unknown turbulence transportMethod: "+method);
        turbulence.stageConservativeTransport=method=="stageConservative";
    }
    value = jsonValueAfterKey(text, "turbulentPrandtl");
    if (!value.empty()) {
        turbulence.turbulentPrandtl =
            valueDoubleValue(value, turbulence.turbulentPrandtl);
    }
    value = jsonValueAfterKey(text, "coupleMomentum");
    if (!value.empty()) {
        turbulence.coupleMomentum =
            valueBoolValue(value, turbulence.coupleMomentum);
    }
    value = jsonValueAfterKey(text, "coupleEnergy");
    if (!value.empty()) {
        turbulence.coupleEnergy =
            valueBoolValue(value, turbulence.coupleEnergy);
    }
}

void CaseAdapter::decodeThermoDynamicsSection() {
    if (!sections_.hasThermoDynamics) return;
    const std::string context = "models/thermoDynamics.yaml";
    const P* text = &sections_.thermoDynamics;

    namespace MP = Physics::Multiphase;
    const P* block = jsonBlockAfterKey(text, "temperature");
    if (!block) block = jsonBlockAfterKey(text, "thermal");
    if (!block) block = jsonBlockAfterKey(text, "thermophysical");
    if (!block) block = text;

    MP::TemperatureOptions temperature =
        multiPhaseLoaded_ ? multiPhaseConfig_.temperature
                          : MP::TemperatureOptions{};

    std::string value = jsonValueAfterKey(block, "enabled");
    if (value.empty()) value = jsonValueAfterKey(block, "enable");
    if (value.empty()) value = jsonValueAfterKey(block, "solve");
    if (value.empty()) value = jsonValueAfterKey(block, "solveTemperature");
    if (value.empty()) value = jsonValueAfterKey(block, "solveEnergy");
    const bool explicitEnabled = !value.empty();
    if (explicitEnabled) {
        temperature.enabled = valueBoolValue(value, temperature.enabled);
    }

    value = jsonValueAfterKey(block, "fieldName");
    if (value.empty()) value = jsonValueAfterKey(block, "name");
    if (!value.empty()) temperature.fieldName = valueUnquote(value);

    value = jsonValueAfterKey(block, "defaultValue");
    if (value.empty()) value = jsonValueAfterKey(block, "defaultTemperature");
    if (value.empty()) value = jsonValueAfterKey(block, "T0");
    if (!value.empty()) {
        temperature.defaultValue =
            valueDoubleValue(value, temperature.defaultValue);
    }

    value = jsonValueAfterKey(block, "transport");
    if (value.empty()) value = jsonValueAfterKey(block, "transportEnabled");
    if (value.empty()) value = jsonValueAfterKey(block, "convection");
    if (value.empty()) value = jsonValueAfterKey(block, "convectionEnabled");
    if (!value.empty()) {
        temperature.transportEnabled =
            valueBoolValue(value, temperature.transportEnabled);
    }

    value = jsonValueAfterKey(block, "diffusion");
    if (value.empty()) value = jsonValueAfterKey(block, "diffusionEnabled");
    if (value.empty()) value = jsonValueAfterKey(block, "conduction");
    if (value.empty()) value = jsonValueAfterKey(block, "heatConduction");
    if (!value.empty()) {
        temperature.diffusionEnabled =
            valueBoolValue(value, temperature.diffusionEnabled);
    }

    value = jsonValueAfterKey(block, "diffusivity");
    if (value.empty()) value = jsonValueAfterKey(block, "thermalDiffusivity");
    if (value.empty()) value = jsonValueAfterKey(block, "alpha");
    if (value.empty()) value = jsonValueAfterKey(block, "DT");
    if (value.empty()) value = jsonValueAfterKey(block, "D");
    if (!value.empty()) {
        temperature.diffusivity =
            valueDoubleValue(value, temperature.diffusivity);
        temperature.diffusionEnabled = temperature.diffusivity > 0.0;
    }

    value = jsonValueAfterKey(block, "bounded");
    if (!value.empty()) {
        temperature.bounded =
            valueBoolValue(value, temperature.bounded);
    }
    value = jsonValueAfterKey(block, "boundMode");
    if (value.empty()) value = jsonValueAfterKey(block, "boundsMode");
    if (!value.empty()) temperature.boundMode = valueUnquote(value);
    value = jsonValueAfterKey(block, "lowerBound");
    if (value.empty()) value = jsonValueAfterKey(block, "min");
    if (!value.empty()) {
        temperature.lowerBound =
            valueDoubleValue(value, temperature.lowerBound);
    }
    value = jsonValueAfterKey(block, "upperBound");
    if (value.empty()) value = jsonValueAfterKey(block, "max");
    if (!value.empty()) {
        temperature.upperBound =
            valueDoubleValue(value, temperature.upperBound);
    }
    value = jsonValueAfterKey(block, "bounds");
    const std::vector<double> bounds = valueNumbers(value);
    if (bounds.size() >= 2) {
        temperature.lowerBound = bounds[0];
        temperature.upperBound = bounds[1];
    }

    if (!explicitEnabled && !temperature.diffusionEnabled) return;
    if (!temperature.enabled) {
        if (multiPhaseLoaded_) multiPhaseConfig_.temperature = temperature;
        return;
    }

    if (!multiPhaseLoaded_) {
        multiPhaseConfig_ = MP::MultiPhaseConfig{};
        multiPhaseConfig_.enabled = true;
        multiPhaseConfig_.type = "thermal";
        multiPhaseLoaded_ = true;
    }
    multiPhaseConfig_.temperature = temperature;
    try {
        MP::validateMultiPhaseConfig(multiPhaseConfig_, context);
    } catch (const std::exception& e) {
        std::cerr << "[SF FATAL] " << e.what() << std::endl;
        std::exit(1);
    }
}


void CaseAdapter::decodeILWSection() {
    if (!sections_.hasIlw) return;
    const std::string context = "models/ILW.yaml";
    broadcast("Loading ILW config: ", context);
    const P* text = &sections_.ilw;
    FDM::BoundaryConfig& boundaries = caseConfig_.solver.boundaries;
    FDM::IBMConfig& ibm = caseConfig_.solver.ibm;
    int ilwOrder = boundaries.ilwOrder;

    {
        const P* ilw = jsonBlockAfterKey(text, "ILW");
        if (!ilw) ilw = text;
        std::string value = jsonValueAfterKey(ilw, "enabled");
        if (value.empty()) value = jsonValueAfterKey(ilw, "enable");
        boundaries.ilwEnabled = !value.empty()
            && valueBoolValue(value, false);
        if (!boundaries.ilwEnabled) {
            ilwOrder = 0;
            boundaries.ilwOrder = 0;
            broadcast("ILW config: ", "disabled");
            return;
        }

        value = jsonValueAfterKey(ilw, "value");
        if (value.empty()) value = jsonValueAfterKey(ilw, "order");
        if (value.empty()) value = jsonValueAfterKey(ilw, "accuracy");
        if (value.empty()) {
            fatalILWConfig(
                "models/ILW.yaml has enabled=true but no value/order.");
        }
        if (!value.empty()) {
            ilwOrder = valueIntValue(value, ilwOrder);
            if (!FDM::isSupportedILWOrder(ilwOrder)) {
                fatalILWConfig("ILW value supports only 0, 3, 5, 7, or 9; got "
                               + std::to_string(ilwOrder) + ".");
            }
        }

        const P* normalSearch = jsonBlockAfterKey(text, "normalSearch");
        if (!normalSearch) normalSearch = jsonBlockAfterKey(ilw, "normalSearch");
        auto parseFoamAngle = [&](const P* block,
                                  const std::string& sectionName) {
            if (!block) return;
            std::string angleValue = jsonValueAfterKey(block, "angle");
            if (angleValue.empty()) angleValue = jsonValueAfterKey(block, "normalAngle");
            if (angleValue.empty()) angleValue = jsonValueAfterKey(block, "normalAngleDegrees");
            if (angleValue.empty()) return;
            const double angle =
                valueDoubleValue(angleValue, ibm.normalAngleDegrees);
            if (!std::isfinite(angle) || angle < 0.0 || angle > 90.0) {
                fatalILWConfig(sectionName
                               + ".angle must be in [0, 90] degrees; got "
                               + std::to_string(angle) + ".");
            }
            ibm.normalAngleDegrees = angle;
        };
        parseFoamAngle(ilw, "ILW");
        parseFoamAngle(normalSearch, "normalSearch");

        value = jsonValueAfterKey(normalSearch, "minLayers");
        if (!value.empty()) {
            ibm.normalSearchMinLayers =
                valueIntValue(value, ibm.normalSearchMinLayers);
        }
        value = jsonValueAfterKey(normalSearch, "maxLayers");
        if (!value.empty()) {
            ibm.normalSearchMaxLayers =
                valueIntValue(value, ibm.normalSearchMaxLayers);
        }
        value = jsonValueAfterKey(normalSearch, "targetCandidates");
        if (!value.empty()) {
            ibm.normalSearchTargetCandidates =
                valueIntValue(value, ibm.normalSearchTargetCandidates);
        }
        value = jsonValueAfterKey(normalSearch, "keepSamples");
        if (!value.empty()) {
            ibm.normalSearchKeepSamples =
                valueIntValue(value, ibm.normalSearchKeepSamples);
        }

        if (ibm.normalSearchMinLayers < 1) {
            fatalILWConfig("normalSearch.minLayers must be >= 1; got "
                           + std::to_string(ibm.normalSearchMinLayers) + ".");
        }
        if (ibm.normalSearchMaxLayers > 0
            && ibm.normalSearchMaxLayers < ibm.normalSearchMinLayers) {
            fatalILWConfig("normalSearch.maxLayers must be >= minLayers or <= 0 for auto; got "
                           + std::to_string(ibm.normalSearchMaxLayers) + ".");
        }
        if (ibm.normalSearchTargetCandidates < 1) {
            fatalILWConfig("normalSearch.targetCandidates must be >= 1; got "
                           + std::to_string(ibm.normalSearchTargetCandidates) + ".");
        }
        if (ibm.normalSearchKeepSamples < 1) {
            fatalILWConfig("normalSearch.keepSamples must be >= 1; got "
                           + std::to_string(ibm.normalSearchKeepSamples) + ".");
        }
        if (ibm.normalSearchKeepSamples > ibm.normalSearchTargetCandidates) {
            fatalILWConfig("normalSearch.keepSamples must be <= targetCandidates; got "
                           + std::to_string(ibm.normalSearchKeepSamples)
                           + " > "
                           + std::to_string(ibm.normalSearchTargetCandidates)
                           + ".");
        }

        boundaries.ilwEnabled = (ilwOrder > 0);
        boundaries.ilwOrder = ilwOrder;
        broadcast("ILW config: ",
                  "value=" + std::to_string(ilwOrder)
                  + ", normalAngle="
                  + std::to_string(ibm.normalAngleDegrees)
                  + " deg");
        return;
    }
}

void CaseAdapter::decodeSourceSettings() {
    auto& sources=caseConfig_.solver.sources;
    sources.gravity.clear();sources.rotating.clear();
    sources.wallHeat=wallHeatBoundarySettings_;
    // Numerical source names do not provide missing model parameters or override model activity.
    const auto requested=sourceScheme_;
    (void)FDM::parseSourceKinds(requested);
    if(sourceSchemeMentions(requested,"gravity") && !sections_.hasGravity)
        throw std::runtime_error("Selected gravity source lacks a registered gravity model.");
    if((sourceSchemeMentions(requested,"mrf") || sourceSchemeMentions(requested,"rotating")) && !sections_.hasMrf)
        throw std::runtime_error("Selected MRF source lacks a registered MRF model.");
    if(sections_.hasGravity)decodeGravitySection();
    if(sections_.hasMrf)decodeMRFSection();
    if(sections_.hasWallHeat)decodeWallHeatSection();
    if(sourceSchemeMentions(requested,"wallHeat") && !sections_.hasWallHeat
        && sources.wallHeat.empty() && multiPhaseConfig_.phaseChange.wallBoiling.empty())
        throw std::runtime_error("Selected wallHeat source lacks a registered model or active wall boundary.");
    sourceScheme_.clear();
    if(!sources.gravity.empty())appendSourceSchemeToken(sourceScheme_,"Gravity");
    if(!sources.rotating.empty())appendSourceSchemeToken(sourceScheme_,"MRF");
    if(!sources.wallHeat.empty())appendSourceSchemeToken(sourceScheme_,"WallHeat");
}

void CaseAdapter::decodeGravitySection() {
    if (!sections_.hasGravity) return;
    const std::string context = "models/gravity.yaml";
    const P* text = &sections_.gravity;
    if(!modelEnabled(text,{"gravity"}))return;

    {
        const P* block = jsonBlockAfterKey(text, "gravity");
        if (!block) block = text;
        const auto acceleration=[&](const P* parameters) {
            if(parameters->contains("value") && parameters->contains("g"))
                throw std::runtime_error("Conflicting explicit gravity value/g selections.");
            if(parameters->contains("value"))return token(parameters->at("value"));
            if(parameters->contains("g"))return token(parameters->at("g"));
            return std::string();
        };
        std::string value=acceleration(block);
        Vector3 g;
        if (!value.empty()) {
            if(!valueVectorValue(value,&g) || !std::isfinite(g[0]) || !std::isfinite(g[1]) || !std::isfinite(g[2]))
                throw std::runtime_error("Invalid gravity acceleration vector.");
            ZoneVectorSetting setting;
            setting.zone = "all";
            setting.value = g;
            caseConfig_.solver.sources.gravity.push_back(setting);
        }
        for (const auto& [zone, zoneBlock] : jsonChildBlocks(block)) {
            value=acceleration(zoneBlock);
            if(!valueVectorValue(value,&g) || !std::isfinite(g[0]) || !std::isfinite(g[1]) || !std::isfinite(g[2]))
                throw std::runtime_error("Invalid gravity acceleration for zone: "+zone);
            ZoneVectorSetting setting;
            setting.zone = zone;
            setting.value = g;
            caseConfig_.solver.sources.gravity.push_back(setting);
        }
        if(caseConfig_.solver.sources.gravity.empty())throw std::runtime_error("Enabled gravity model requires a finite value/g vector or zone contributions.");
        return;
    }
}

void CaseAdapter::decodeMRFSection() {
    if (!sections_.hasMrf) return;
    const std::string context = "models/MRF.yaml";
    const P* text = &sections_.mrf;
    if(!modelEnabled(text,{"MRF","rotate"}))return;

    {
        auto parseRotateBlock = [&](const P* block,
                                    const std::string& fallbackZone) {
            if (!block) return;
            RotatingSetting rotate;
            rotate.zone = fallbackZone.empty() ? "all" : fallbackZone;
            std::string value = jsonValueAfterKey(block, "zone");
            if (value.empty()) value = jsonValueAfterKey(block, "set");
            if (!value.empty()) rotate.zone = valueUnquote(value);
            Vector3 vec;
            value = jsonValueAfterKey(block, "center");
            if (valueVectorValue(value, &vec)) rotate.center = vec;
            value = jsonValueAfterKey(block, "axis");
            if (valueVectorValue(value, &vec)) rotate.axis = normalize(vec);
            value = jsonValueAfterKey(block, "omega");
            if (value.empty()) value = jsonValueAfterKey(block, "speed");
            if (!value.empty()) rotate.omega = valueDoubleValue(value, rotate.omega);
            value = jsonValueAfterKey(block, "velocity");
            if (valueVectorValue(value, &vec)) {
                rotate.velocity = vec;
                rotate.hasVelocity = true;
            }
            caseConfig_.solver.sources.rotating.push_back(rotate);
        };

        const P* block = jsonBlockAfterKey(text, "MRF");
        if (!block) block = jsonBlockAfterKey(text, "rotate");
        if (block) parseRotateBlock(block, "all");
        const P* zones = jsonBlockAfterKey(text, "zones");
        for (const auto& [zone, zoneBlock] : jsonChildBlocks(zones)) {
            parseRotateBlock(zoneBlock, zone);
        }
        if(caseConfig_.solver.sources.rotating.empty())throw std::runtime_error("Enabled MRF model requires MRF/rotate or zones parameters.");
        return;
    }
}

void CaseAdapter::decodeWallHeatSection() {
    if (!sections_.hasWallHeat) return;
    const std::string context = "models/wallHeat.yaml";
    const P* text = &sections_.wallHeat;
    if(!modelEnabled(text,{"wallHeatSource","wallHeat"}))return;

    const P* block = jsonBlockAfterKey(text, "wallHeatSource");
    if (!block) block = jsonBlockAfterKey(text, "wallHeat");
    if (!block) block = text;

    std::string value = jsonValueAfterKey(block, "enabled");
    if (value.empty()) value = jsonValueAfterKey(block, "enable");
    const bool enabled = valueBoolValue(value, true);
    if (!enabled) return;

    std::vector<WallHeatSetting>& wallHeatSettings =
        caseConfig_.solver.sources.wallHeat;
    const size_t before = wallHeatSettings.size();

    auto patchNames = [](const std::string& text,
                         const std::string& fallback) {
        std::vector<std::string> patches;
        if (!text.empty()) patches = valueWordList(text);
        if (patches.empty() && !fallback.empty()) patches.push_back(fallback);
        return patches;
    };

    auto parseHeatBlock = [&](const P* heatBlock,
                              const std::string& fallbackPatch) {
        if (!heatBlock) return;

        std::string local = jsonValueAfterKey(heatBlock, "enabled");
        if (local.empty()) local = jsonValueAfterKey(heatBlock, "enable");
        if (local.empty()) local = jsonValueAfterKey(heatBlock, "active");
        if (!valueBoolValue(local, true)) return;

        WallHeatSetting base;
        local = jsonValueAfterKey(heatBlock, "type");
        if (local.empty()) local = jsonValueAfterKey(heatBlock, "mode");
        if (!local.empty()) base.mode = valueUnquote(local);

        local = jsonValueAfterKey(heatBlock, "coupleEnergy");
        if (local.empty()) local = jsonValueAfterKey(heatBlock, "energyCoupling");
        if (!local.empty()) base.coupleEnergy = valueBoolValue(local, true);

        std::string heatFlux = jsonValueAfterKey(heatBlock, "q");
        if (heatFlux.empty()) heatFlux = jsonValueAfterKey(heatBlock, "heatFlux");
        if (heatFlux.empty()) heatFlux = jsonValueAfterKey(heatBlock, "value");
        if (heatFlux.empty()) {
            fatalCaseConfig("WallHeat block for patch '" + fallbackPatch
                            + "' must provide q or heatFlux.");
        }
        base.heatFlux = valueDoubleValue(heatFlux, base.heatFlux);

        std::string patchText = jsonValueAfterKey(heatBlock, "patch");
        if (patchText.empty()) patchText = jsonValueAfterKey(heatBlock, "patches");
        if (patchText.empty()) patchText = jsonValueAfterKey(heatBlock, "set");
        if (patchText.empty()) patchText = jsonValueAfterKey(heatBlock, "zone");
        const std::vector<std::string> patches =
            patchNames(patchText, fallbackPatch);
        if (patches.empty()) {
            fatalCaseConfig("WallHeat block must name patch/patches.");
        }

        for (const std::string& patch : patches) {
            WallHeatSetting setting = base;
            setting.patch = patch;
            wallHeatSettings.push_back(setting);
        }
    };

    if (!jsonValueAfterKey(block, "q").empty()
        || !jsonValueAfterKey(block, "heatFlux").empty()
        || !jsonValueAfterKey(block, "value").empty()) {
        parseHeatBlock(block, "");
    }

    const P* sources = jsonBlockAfterKey(block, "sources");
    for (const auto& [name, heatBlock] : jsonChildBlocks(sources)) {
        parseHeatBlock(heatBlock, name);
    }

    const P* patchesBlock = jsonBlockAfterKey(block, "patches");
    for (const auto& [patch, heatBlock] : jsonChildBlocks(patchesBlock)) {
        parseHeatBlock(heatBlock, patch);
    }

    const P* namedSources = jsonBlockAfterKey(text, "wallHeatSources");
    for (const auto& [name, heatBlock] : jsonChildBlocks(namedSources)) {
        parseHeatBlock(heatBlock, name);
    }

    if (wallHeatSettings.size() == before) {
        fatalCaseConfig("WallHeat source file is enabled but no source "
                        "settings were parsed: " + context);
    }
}

// ============================================================
//  SFM 网格 block 读取在 src/infrastructure/mesh/MultiBlockMesh/ 中实现:
//  每个 #Information 是一个 block 段，同段重复 #tag 会归并。
// ============================================================

// ============================================================
//  VTK 输出 — 重构版: 支持 PVD/VTM/VTS 层级输出
// ============================================================

// ============================================================
//  工具
// ============================================================


namespace {

BCType explicitPhiType(const std::string& rawType,
                       const std::string& fileName,
                       const std::string& context) {

    const std::string key = valueLower(valueUnquote(rawType));
    if (key == "fixedvalue" || key == "fixed_value") return FIXED_VALUE;
    if (key == "zerogradient" || key == "zero_gradient") return ZERO_GRADIENT;
    if (key == "empty") return EMPTY;
    fatalCaseConfig(
        fileName + ": " + context
        + " must explicitly use fixedValue, zeroGradient, or empty.");
    return FIXED_VALUE;
}

} // namespace

void CaseAdapter::decodePhiField(const Model::FieldDescriptor& field) {

    namespace MP = Physics::Multiphase;

    multiPhaseConfig_.phiInitialConditions.clear();
    multiPhaseConfig_.phiPlaneInitializers.clear();
    multiPhaseConfig_.phiSphereInitializers.clear();
    multiPhaseConfig_.phiBoundaryConditions.clear();

    double value = std::numeric_limits<double>::quiet_NaN();
    std::string raw = token(field.initial);
    if (raw.empty() || !valueScalarValue(raw, &value)) {
        fatalCaseConfig(
            field.name + ": internalField must explicitly provide a scalar phi value.");
    }
    MP::LevelSetScalarCondition internal;
    internal.setName = "all";
    internal.type = FIXED_VALUE;
    internal.typeDeclared = true;
    internal.value = value;
    multiPhaseConfig_.phiInitialConditions.push_back(internal);

    const P* sets = &field.sets;
    for (const auto& [name, block] : jsonChildBlocks(sets)) {
        const std::string type =
            valueLower(valueUnquote(jsonValueAfterKey(block, "type")));
        if (type == "signeddistanceplane" || type == "signed_distance_plane"
            || type == "planedistance" || type == "plane_distance") {
            MP::LevelSetPlaneInitializer plane;
            plane.setName = name;
            if (!valueVectorValue(jsonValueAfterKey(block, "point"), &plane.point)
                || !valueVectorValue(jsonValueAfterKey(block, "normal"),
                                    &plane.normal)) {
                fatalCaseConfig(
                    field.name + ": signed-distance plane '" + name
                    + "' requires explicit point and normal vectors.");
            }
            raw = jsonValueAfterKey(block, "scale");
            if (raw.empty()) {
                fatalCaseConfig(
                    field.name + ": signed-distance plane '" + name
                    + "' requires an explicit positive scale.");
            }
            plane.scale = valueDoubleValue(raw, plane.scale);
            multiPhaseConfig_.phiPlaneInitializers.push_back(plane);
            continue;
        }
        if (type == "signeddistancesphere" || type == "signed_distance_sphere"
            || type == "spheredistance" || type == "sphere_distance") {
            MP::LevelSetSphereInitializer sphere;
            sphere.setName = name;
            if (!valueVectorValue(jsonValueAfterKey(block, "center"),
                                 &sphere.center)) {
                fatalCaseConfig(
                    field.name + ": signed-distance sphere '" + name
                    + "' requires an explicit center vector.");
            }
            raw = jsonValueAfterKey(block, "radius");
            if (raw.empty()) {
                fatalCaseConfig(
                    field.name + ": signed-distance sphere '" + name
                    + "' requires an explicit positive radius.");
            }
            sphere.radius = valueDoubleValue(raw, sphere.radius);
            raw = jsonValueAfterKey(block, "scale");
            if (raw.empty()) {
                fatalCaseConfig(
                    field.name + ": signed-distance sphere '" + name
                    + "' requires an explicit positive scale.");
            }
            sphere.scale = valueDoubleValue(raw, sphere.scale);
            sphere.insidePhase =
                valueUnquote(jsonValueAfterKey(block, "insidePhase"));
            if (sphere.insidePhase.empty()) {
                fatalCaseConfig(
                    field.name + ": signed-distance sphere '" + name
                    + "' requires an explicit insidePhase.");
            }
            multiPhaseConfig_.phiSphereInitializers.push_back(sphere);
            continue;
        }

        raw = jsonValueAfterKey(block, "value");
        if (!valueScalarValue(raw, &value)) {
            fatalCaseConfig(
                field.name + ": internalSets entry '" + name
                + "' requires an explicit scalar value.");
        }
        MP::LevelSetScalarCondition condition;
        condition.setName = name;
        condition.type = explicitPhiType(
            jsonValueAfterKey(block, "type"), field.name,
            "internalSets entry '" + name + "'");
        condition.typeDeclared = true;
        condition.value = value;
        multiPhaseConfig_.phiInitialConditions.push_back(condition);
    }

    const P* patches = &field.boundaries;
    for (const auto& [name, block] : jsonChildBlocks(patches)) {
        MP::LevelSetScalarCondition condition;
        condition.setName = name;
        condition.type = explicitPhiType(
            jsonValueAfterKey(block, "type"), field.name,
            "boundaryField entry '" + name + "'");
        condition.typeDeclared = true;
        if (condition.type == FIXED_VALUE) {
            raw = jsonValueAfterKey(block, "value");
            if (!valueScalarValue(raw, &value)) {
                fatalCaseConfig(
                    field.name + ": fixedValue boundary '" + name
                    + "' requires an explicit scalar value.");
            }
            condition.value = value;
        }
        multiPhaseConfig_.phiBoundaryConditions.push_back(condition);
    }
}


namespace {

struct IBMGeometryToken {
    std::string value;
    bool explicitPath = false;
};

/// @brief 解析 IBM STL 列表并保留“是否带引号”的路径语义。
std::vector<IBMGeometryToken> ibmGeometryTokens(const std::string& text) {
    std::vector<IBMGeometryToken> tokens;
    std::size_t pos = 0;
    auto separator = [](char c) {
        return std::isspace(static_cast<unsigned char>(c))
            || c == '(' || c == ')' || c == '[' || c == ']'
            || c == ',' || c == ';';
    };

    while (pos < text.size()) {
        while (pos < text.size() && separator(text[pos])) ++pos;
        if (pos >= text.size()) break;

        IBMGeometryToken token;
        if (text[pos] == '"' || text[pos] == '\'') {
            token.explicitPath = true;
            const char quote = text[pos++];
            const std::size_t begin = pos;
            while (pos < text.size() && text[pos] != quote) ++pos;
            if (pos >= text.size()) {
                fatalIBMConfig(
                    "IBM geometryFiles contains an unterminated quoted path.");
            }
            token.value = text.substr(begin, pos - begin);
            ++pos;
        } else {
            const std::size_t begin = pos;
            while (pos < text.size() && !separator(text[pos])) ++pos;
            token.value = text.substr(begin, pos - begin);
        }
        if (!token.value.empty()) tokens.push_back(std::move(token));
    }
    return tokens;
}

} // namespace

void CaseAdapter::addIBMGeometryFile(const std::string& fileName,
                            bool explicitPath) {
    std::string clean = stripTokenQuotes(fileName);
    if (clean.empty()) return;

    // geometry 路径是 case-relative 语义，不按历史目录改写；未给出路径即 fail fast。
    if (!explicitPath && !isPlainRelativeFileName(clean)) {
        fatalIBMConfig(
            "IBM geometry entry must be a case-relative path such as "
            "\"mesh/geometry/cylinder.stl\". Got: " + clean);
    }

    const std::string absPath = resolveCaseFile(caseDir_, clean);
    if (!fileExists(absPath)) {
        fatalIBMConfig("IBM geometry file not found: " + absPath);
    }

    if (std::find(ibmGeometryFiles_.begin(), ibmGeometryFiles_.end(), clean)
        == ibmGeometryFiles_.end()) {
        ibmGeometryFiles_.push_back(clean);
    }
}

void CaseAdapter::decodeIBMSection() {
    if (!sections_.hasIbm) return;
    const std::string context = "models/IBM.yaml";
    const P* text = &sections_.ibm;

    const P* block = jsonBlockAfterKey(text, "IBM");
    if (!block) block = text;
    std::string value = jsonValueAfterKey(block, "enabled");
    if (value.empty()) value = jsonValueAfterKey(block, "enable");
    caseConfig_.solver.ibm.enabled = !value.empty()
        && valueBoolValue(value, false);
    if (!caseConfig_.solver.ibm.enabled) {
        ibmGeometryFiles_.clear();
        return;
    }

    value = jsonValueAfterKey(block, "type");
    try {
        ibmMethod_ = value.empty()
            ? FDM::IBMMethod::Ghost
            : FDM::parseIBMMethod(valueUnquote(value));
    } catch (const std::exception& error) {
        fatalIBMConfig(error.what());
    }
    const auto wall=valueUnquote(jsonValueAfterKey(block,"wallClosure"));
    if (wall.empty() || wall=="eulerSlip")
        caseConfig_.solver.ibm.wallClosure=FDM::ImmersedWallClosure::EulerSlip;
    else if (wall=="stationaryNoSlipAdiabatic")
        caseConfig_.solver.ibm.wallClosure=FDM::ImmersedWallClosure::StationaryNoSlipAdiabatic;
    else fatalIBMConfig("Unknown immersed wallClosure: "+wall);
    if (ibmMethod_!=FDM::IBMMethod::Ghost && !wall.empty())
        fatalIBMConfig("wallClosure belongs to Ghost boundary; forcing does not implement this wall closure.");
    if (ibmMethod_ == FDM::IBMMethod::VariationalForcing) {
        const P* forcing = jsonBlockAfterKey(block, "forcing");
        if (!forcing) {
            fatalIBMConfig(
                "variationalForcing requires an explicit IBM/forcing block.");
        }
        auto required = [&](const std::string& key) {
            const std::string raw = jsonValueAfterKey(forcing, key);
            if (raw.empty()) {
                fatalIBMConfig(
                    "variationalForcing requires explicit forcing."
                    + key + ".");
            }
            return raw;
        };
        const std::string formulation =
            FDM::normalizeToken(required("formulation"));
        if (formulation != "variationaldlm") {
            fatalIBMConfig(
                "forcing.formulation must be variationalDLM; got '"
                + formulation + "'.");
        }
        const std::string algorithmValue = required("algorithm");
        try {
            ibmForcingConfig_.algorithm =
                FDM::parseIBMForcingAlgorithm(algorithmValue);
        } catch (const std::exception& error) {
            fatalIBMConfig(error.what());
        }
        try {
            const std::string domain =
                jsonValueAfterKey(forcing, "constraintDomain");
            const std::string support =
                jsonValueAfterKey(forcing, "constraintSupport");
            if (domain.empty() && support.empty()) {
                fatalIBMConfig(
                    "variationalForcing requires explicit forcing."
                    "constraintSupport (or legacy constraintDomain).");
            }
            if (!domain.empty()) {
                ibmForcingConfig_.constraintDomain =
                    FDM::parseIBMConstraintDomain(domain);
            }
            if (!support.empty()) {
                ibmForcingConfig_.constraintSupport =
                    FDM::parseIBMConstraintSupport(support);
            } else {
                ibmForcingConfig_.constraintSupport =
                    ibmForcingConfig_.constraintDomain
                        == FDM::IBMConstraintDomain::Surface
                    ? FDM::IBMConstraintSupport::Surface
                    : FDM::IBMConstraintSupport::Body;
            }
            // 新版 support 仍需填充 legacy 字段，因为老的 forcing 实现和
            // 外部观察器可能读取 constraintDomain；这里做一次显式同步，
            // 不让默认的 volume 值悄悄覆盖用户选择的 surface。
            if (domain.empty()
                && ibmForcingConfig_.constraintSupport
                    != FDM::IBMConstraintSupport::SurfaceAndBody) {
                ibmForcingConfig_.constraintDomain =
                    ibmForcingConfig_.constraintSupport
                        == FDM::IBMConstraintSupport::Surface
                    ? FDM::IBMConstraintDomain::Surface
                    : FDM::IBMConstraintDomain::Volume;
            }
            if (!domain.empty() && !support.empty()) {
                const FDM::IBMConstraintSupport domainSupport =
                    ibmForcingConfig_.constraintDomain
                        == FDM::IBMConstraintDomain::Surface
                    ? FDM::IBMConstraintSupport::Surface
                    : FDM::IBMConstraintSupport::Body;
                if (ibmForcingConfig_.constraintSupport != domainSupport
                    && ibmForcingConfig_.constraintSupport
                        != FDM::IBMConstraintSupport::SurfaceAndBody) {
                    fatalIBMConfig(
                        "forcing.constraintDomain and constraintSupport "
                        "select different domains.");
                }
            }
            const std::string coupling =
                jsonValueAfterKey(forcing, "coupling");
            const std::string enforcement =
                jsonValueAfterKey(forcing, "enforcement");
            if (coupling.empty() && enforcement.empty()) {
                fatalIBMConfig(
                    "forcing requires explicit coupling or enforcement.");
            }
            if (!enforcement.empty()) {
                ibmForcingConfig_.enforcement =
                    FDM::parseIBMEnforcement(enforcement);
            }
            if (!coupling.empty()) {
                ibmForcingConfig_.coupling =
                    FDM::parseIBMConstraintCoupling(coupling);
            } else {
                ibmForcingConfig_.coupling =
                    ibmForcingConfig_.enforcement
                        == FDM::IBMEnforcement::MonolithicKKT
                    ? FDM::IBMConstraintCoupling::MonolithicKKT
                    : FDM::IBMConstraintCoupling::IncrementalProjection;
            }
            if (enforcement.empty()) {
                ibmForcingConfig_.enforcement =
                    ibmForcingConfig_.coupling
                        == FDM::IBMConstraintCoupling::MonolithicKKT
                    ? FDM::IBMEnforcement::MonolithicKKT
                    : FDM::IBMEnforcement::FractionalDLM;
            }
            if ((!coupling.empty() && !enforcement.empty())
                && ((ibmForcingConfig_.coupling
                         == FDM::IBMConstraintCoupling::MonolithicKKT)
                    != (ibmForcingConfig_.enforcement
                            == FDM::IBMEnforcement::MonolithicKKT))) {
                fatalIBMConfig(
                    "forcing coupling and enforcement select different "
                    "time-level strategies.");
            }
            const std::string representation =
                jsonValueAfterKey(forcing, "representation");
            if (!representation.empty()) {
                ibmForcingConfig_.representation =
                    FDM::parseIBMRepresentation(representation);
            } else if (ibmForcingConfig_.constraintSupport
                       == FDM::IBMConstraintSupport::Surface) {
                ibmForcingConfig_.representation =
                    FDM::IBMRepresentation::DiffuseKernel;
            } else if (ibmForcingConfig_.constraintSupport
                       == FDM::IBMConstraintSupport::Body) {
                ibmForcingConfig_.representation =
                    FDM::IBMRepresentation::EulerianMask;
            } else {
                fatalIBMConfig(
                    "constraintSupport surfaceAndBody requires explicit "
                    "forcing.representation.");
            }
            const std::string motion = jsonValueAfterKey(forcing, "motion");
            const std::string solidModel =
                jsonValueAfterKey(forcing, "solidModel");
            if (motion.empty() && solidModel.empty()) {
                fatalIBMConfig(
                    "forcing requires explicit motion or solidModel.");
            }
            if (!solidModel.empty()) {
                ibmForcingConfig_.solidModel =
                    FDM::parseIBMSolidModel(solidModel);
            }
            if (!motion.empty()) {
                ibmForcingConfig_.motion =
                    FDM::parseIBMSolidMotion(motion);
            } else {
                ibmForcingConfig_.motion =
                    ibmForcingConfig_.solidModel
                        == FDM::IBMSolidModel::CoupledRigid
                    || ibmForcingConfig_.solidModel
                        == FDM::IBMSolidModel::SelfPropelledRigid
                    ? FDM::IBMSolidMotion::CoupledRigid
                    : FDM::IBMSolidMotion::PrescribedRigid;
            }
            if (solidModel.empty()) {
                ibmForcingConfig_.solidModel =
                    ibmForcingConfig_.motion
                        == FDM::IBMSolidMotion::CoupledRigid
                    ? FDM::IBMSolidModel::CoupledRigid
                    : FDM::IBMSolidModel::Prescribed;
            }
            ibmForcingConfig_.energyCoupling =
                FDM::parseIBMEnergyCoupling(required("energyCoupling"));
        } catch (const std::exception& error) {
            fatalIBMConfig(error.what());
        }
        if (!valueVectorValue(required("centerOfMass"),
                             &ibmForcingConfig_.centerOfMass)
            || !valueVectorValue(required("linearVelocity"),
                                &ibmForcingConfig_.linearVelocity)
            || !valueVectorValue(required("angularVelocity"),
                                &ibmForcingConfig_.angularVelocity)) {
            fatalIBMConfig(
                "forcing centerOfMass/linearVelocity/angularVelocity must "
                "each be a three-component SI vector.");
        }
        try {
            ibmForcingConfig_.rigidMotionMode =
                FDM::parseIBMRigidMotionMode(required("motionMode"));
        } catch (const std::exception& error) {
            fatalIBMConfig(error.what());
        }
        ibmForcingConfig_.fluidPorts.clear();
        std::string ports = jsonValueAfterKey(forcing, "fluidPorts");
        if (ports.empty()) ports = jsonValueAfterKey(forcing, "phases");
        for (const std::string& rawPort : FDM::splitList(ports)) {
            const std::string port = stripTokenQuotes(rawPort);
            if (!port.empty()) ibmForcingConfig_.fluidPorts.push_back(port);
        }
        ibmForcingConfig_.constraintTolerance = valueDoubleValue(
            required("constraintTolerance"),
            ibmForcingConfig_.constraintTolerance);
        if (!std::isfinite(ibmForcingConfig_.constraintTolerance)
            || ibmForcingConfig_.constraintTolerance <= 0.0) {
            fatalIBMConfig(
                "forcing.constraintTolerance must be finite and > 0 m/s.");
        }
        if (ibmForcingConfig_.algorithm
            == FDM::IBMForcingAlgorithm::VelocityForcingFTS) {
            const P* solver = jsonBlockAfterKey(forcing, "constraintSolver");
            if (!solver) {
                fatalIBMConfig(
                    "velocityForcingFTS requires explicit "
                    "forcing.constraintSolver controls.");
            }
            const std::string type = valueUnquote(
                jsonValueAfterKey(solver,"type"));
            const std::string maxIterations =
                jsonValueAfterKey(solver,"maxIterations");
            const std::string relativeTolerance =
                jsonValueAfterKey(solver,"relativeTolerance");
            if (FDM::normalizeToken(type) != "matrixfreecg"
                || maxIterations.empty() || relativeTolerance.empty()) {
                fatalIBMConfig(
                    "constraintSolver requires type matrixFreeCG, explicit "
                    "maxIterations and relativeTolerance.");
            }
            ibmForcingConfig_.constraintSolverMaxIterations =
                valueIntValue(maxIterations,0);
            ibmForcingConfig_.constraintSolverRelativeTolerance =
                valueDoubleValue(relativeTolerance,
                    ibmForcingConfig_.constraintSolverRelativeTolerance);
        }
        if (ibmForcingConfig_.constraintSupport
            == FDM::IBMConstraintSupport::Surface
            || ibmForcingConfig_.constraintSupport
                == FDM::IBMConstraintSupport::SurfaceAndBody) {
            const P* surface = jsonBlockAfterKey(forcing, "surfaceOperator");
            if (!surface) {
                fatalIBMConfig(
                    "surface constraint requires forcing.surfaceOperator.");
            }
            auto surfaceRequired = [&](const std::string& key) {
                const std::string raw = jsonValueAfterKey(surface, key);
                if (raw.empty()) {
                    fatalIBMConfig(
                        "surfaceOperator requires explicit " + key + ".");
                }
                return valueUnquote(raw);
            };
            ibmForcingConfig_.surfaceKernel =
                surfaceRequired("kernel");
            ibmForcingConfig_.surfaceQuadrature =
                surfaceRequired("quadrature");
            const auto normalization=FDM::normalizeToken(surfaceRequired("normalization"));
            if (normalization=="partitionofunity")
                ibmForcingConfig_.surfaceNormalization=FDM::IBMSurfaceNormalization::PartitionOfUnity;
            else if (normalization=="linearreproducing")
                ibmForcingConfig_.surfaceNormalization=FDM::IBMSurfaceNormalization::LinearReproducing;
            else fatalIBMConfig("surfaceOperator normalization must be partitionOfUnity or linearReproducing.");
            ibmForcingConfig_.surfaceSpreading =
                surfaceRequired("spreading");
            ibmForcingConfig_.surfaceSupportRadius = valueDoubleValue(
                surfaceRequired("supportRadius"),
                ibmForcingConfig_.surfaceSupportRadius);
            if (FDM::normalizeToken(ibmForcingConfig_.surfaceKernel)
                    != "wendlandc2"
                || FDM::normalizeToken(
                       ibmForcingConfig_.surfaceQuadrature)
                    != "trianglecentroid"
                || FDM::normalizeToken(
                       ibmForcingConfig_.surfaceSpreading) != "adjoint"
                || !std::isfinite(
                    ibmForcingConfig_.surfaceSupportRadius)
                || ibmForcingConfig_.surfaceSupportRadius <= 0.0) {
                fatalIBMConfig(
                    "surfaceOperator requires kernel wendlandC2, quadrature "
                    "triangleCentroid, normalization partitionOfUnity or linearReproducing, "
                    "spreading adjoint, and supportRadius > 0 m.");
            }
        }
        const std::string penalty =
            jsonValueAfterKey(forcing, "penaltyCoefficient");
        if (ibmForcingConfig_.enforcement == FDM::IBMEnforcement::BrinkmanPenalty) {
            if (penalty.empty()) {
                fatalIBMConfig(
                    "brinkmanPenalty requires explicit "
                    "forcing.penaltyCoefficient.");
            }
            ibmForcingConfig_.penaltyCoefficient = valueDoubleValue(
                penalty, ibmForcingConfig_.penaltyCoefficient);
        }
        const std::string augmentation =
            jsonValueAfterKey(forcing, "augmentationCoefficient");
        if (ibmForcingConfig_.algorithm
            == FDM::IBMForcingAlgorithm::DFMAugmentedLagrangian) {
            if (augmentation.empty()) {
                fatalIBMConfig(
                    "dfmAugmentedLagrangian requires explicit "
                    "forcing.augmentationCoefficient.");
            }
            ibmForcingConfig_.augmentationCoefficient=valueDoubleValue(
                augmentation,
                ibmForcingConfig_.augmentationCoefficient);
        }
        if (ibmForcingConfig_.motion == FDM::IBMSolidMotion::CoupledRigid) {
            const P* rigid = jsonBlockAfterKey(forcing, "rigidBody");
            if (!rigid) {
                fatalIBMConfig(
                    "coupledRigid requires forcing.rigidBody.");
            }
            auto rigidRequired = [&](const std::string& key) {
                const std::string raw = jsonValueAfterKey(rigid, key);
                if (raw.empty()) {
                    fatalIBMConfig(
                        "rigidBody requires explicit " + key + ".");
                }
                return raw;
            };
            ibmForcingConfig_.rigidMass = valueDoubleValue(
                rigidRequired("mass"), ibmForcingConfig_.rigidMass);
            if (!valueVectorValue(rigidRequired("momentOfInertia"),
                                 &ibmForcingConfig_.principalInertia)
                || !valueVectorValue(rigidRequired("externalForce"),
                                    &ibmForcingConfig_.externalForce)
                || !valueVectorValue(rigidRequired("externalTorque"),
                                    &ibmForcingConfig_.externalTorque)
                || !std::isfinite(ibmForcingConfig_.rigidMass)
                || ibmForcingConfig_.rigidMass <= 0.0
                || !std::isfinite(ibmForcingConfig_.principalInertia.x)
                || !std::isfinite(ibmForcingConfig_.principalInertia.y)
                || !std::isfinite(ibmForcingConfig_.principalInertia.z)
                || ibmForcingConfig_.principalInertia.x <= 0.0
                || ibmForcingConfig_.principalInertia.y <= 0.0
                || ibmForcingConfig_.principalInertia.z <= 0.0
                || !std::isfinite(ibmForcingConfig_.externalForce.x)
                || !std::isfinite(ibmForcingConfig_.externalForce.y)
                || !std::isfinite(ibmForcingConfig_.externalForce.z)
                || !std::isfinite(ibmForcingConfig_.externalTorque.x)
                || !std::isfinite(ibmForcingConfig_.externalTorque.y)
                || !std::isfinite(ibmForcingConfig_.externalTorque.z)) {
                fatalIBMConfig(
                    "rigidBody mass must be > 0 kg; momentOfInertia must "
                    "contain three positive kg m2 values; externalForce "
                    "and externalTorque must be SI vectors.");
            }
        }
        try {
            FDM::validateIBMForcingSelection(ibmForcingConfig_);
        } catch (const std::exception& error) {
            fatalIBMConfig(error.what());
        }
    }

    value = jsonValueAfterKey(block, "geometry");
    if (value.empty()) value = jsonValueAfterKey(block, "geometryFiles");
    if (value.empty()) value = jsonValueAfterKey(block, "sheet");
    for (const IBMGeometryToken& file : ibmGeometryTokens(value)) {
        addIBMGeometryFile(file.value, file.explicitPath);
    }

    const P* boundary = jsonBlockAfterKey(block, "boundary");
    value = jsonValueAfterKey(boundary, "geometry");
    if (value.empty()) value = jsonValueAfterKey(boundary, "geometryFiles");
    if (value.empty()) value = jsonValueAfterKey(boundary, "sheet");
    for (const IBMGeometryToken& file : ibmGeometryTokens(value)) {
        addIBMGeometryFile(file.value, file.explicitPath);
    }
    if (ibmGeometryFiles_.empty()) {
        fatalIBMConfig(
            "models/IBM.yaml has enabled=true but no geometryFiles.");
    }
}


void CaseAdapter::decodePhaseChangeSection() {
    if (!sections_.hasPhaseChange) return;
    const std::string context = "models/phaseChange.yaml";
    const P* text = &sections_.phaseChange;

    if (!multiPhaseLoaded_) {
        fatalCaseConfig("models/phaseChange.yaml requires models/multiPhase.yaml "
                        "with a mixture multiPhase configuration: " + context);
    }

    namespace MP = Physics::Multiphase;
    MP::PhaseChangeOptions config = multiPhaseConfig_.phaseChange;

    const P* block = jsonBlockAfterKey(text, "phaseChange");
    if (!block) block = text;
    config.enabled = true;

    auto boolValue = [](const std::string& raw, bool fallback) {
        const std::string value = valueLower(valueUnquote(raw));
        if (value == "ture") return true;
        return valueBoolValue(raw, fallback);
    };

    std::string v = jsonValueAfterKey(block, "enabled");
    if (v.empty()) v = jsonValueAfterKey(block, "phaseChange");
    if (v.empty()) v = jsonValueAfterKey(block, "enable");
    if (!v.empty()) config.enabled = boolValue(v, config.enabled);

    v = jsonValueAfterKey(block, "phaseChangeModel");
    if (v.empty()) {
        v = jsonValueAfterKey(block, "phaseChangeTwoPhaseMixture");
    }
    if (v.empty()) v = jsonValueAfterKey(block, "model");
    if (v.empty()) v = jsonValueAfterKey(block, "type");
    if (!v.empty()) {
        config.model = valueUnquote(v);
    }
    v = jsonValueAfterKey(block, "phases");
    if (!v.empty()) config.phaseNames = valueWordList(v);

    auto modelSubBlock = [&]() -> const P* {
        const std::string normalized =
            Physics::PhaseChange::normalizeModel(config.model);
        const std::string modelKey = sourceTokenKey(config.model);
        for (const auto& [name, child] : jsonChildBlocks(block)) {
            if (Physics::PhaseChange::normalizeModel(name) == normalized
                || sourceTokenKey(name) == modelKey) {
                return child;
            }
        }
        if (normalized == "rpi" || normalized == "wallboiling") {
            const P* child = jsonBlockAfterKey(block, "RPI");
            if (!child) child = jsonBlockAfterKey(block, "wallBoilingModel");
            return child;
        }
        return jsonBlockAfterKey(block, config.model);
    }();

    auto readScalar = [&](const std::string& key,
                          double& target,
                          const std::string& coefficientName = "") {
        std::string value = jsonValueAfterKey(block, key);
        if (value.empty()) return;
        target = valueDoubleValue(value, target);
        if (!coefficientName.empty()) {
            config.coefficients[coefficientName] = target;
        }
    };

    const P* saturation = jsonBlockAfterKey(block, "saturation");
    auto readSaturationScalar = [&](const std::string& key,
                                    double& target) {
        std::string value = jsonValueAfterKey(saturation, key);
        if (value.empty()) value = jsonValueAfterKey(block, key);
        if (!value.empty()) target = valueDoubleValue(value, target);
    };
    readSaturationScalar(
        "saturationTemperature", config.saturationTemperature);
    readSaturationScalar("Tsat", config.saturationTemperature);
    readSaturationScalar(
        "saturationPressure", config.saturationPressure);
    readSaturationScalar("pSat", config.saturationPressure);
    readSaturationScalar("psat", config.saturationPressure);
    readSaturationScalar("latentHeat", config.latentHeat);
    readSaturationScalar("L", config.latentHeat);
    readScalar("evaporationCoefficient", config.evaporationCoefficient);
    readScalar("Ce", config.evaporationCoefficient);
    readScalar("condensationCoefficient", config.condensationCoefficient);
    readScalar("Cc", config.condensationCoefficient);
    config.coefficients["saturationtemperature"] = config.saturationTemperature;
    config.coefficients["saturationpressure"] = config.saturationPressure;
    config.coefficients["latentheat"] = config.latentHeat;
    config.coefficients["evaporationcoefficient"] = config.evaporationCoefficient;
    config.coefficients["condensationcoefficient"] = config.condensationCoefficient;

    v = jsonValueAfterKey(block, "energyCoupling");
    if (!v.empty()) config.energyCoupling = boolValue(v, config.energyCoupling);

    const std::vector<std::string> coefficientKeys = {
        "interfaceAreaDensity",
        "liquidTemperatureGradient",
        "vaporTemperatureGradient",
        "liquidThermalConductivity",
        "vaporThermalConductivity",
        "accommodationCoefficient",
        "molarMass",
        "gasConstant",
        "pSatCoefficient",
        "meltingTemperature",
        "mushyTemperatureWidth",
        "relaxationCoefficient",
        "evaporativeFraction",
        "bubbleInfluenceAreaFraction",
        "departureDiameter",
        "departureFrequency",
        "nucleationSiteDensity",
        "referenceDepartureDiameter",
        "referenceSubcooling",
        "unalA",
        "nearWallLiquidVelocity",
        "unalReferenceVelocity",
        "gravityMagnitude",
        "lemmertChawlaM",
        "lemmertChawlaExponent",
        "averageCavityDensity",
        "contactAngleScale",
        "cavityLengthScale",
        "hibikiDensityFunction",
        "criticalCavityRadius",
        "convectiveHeatTransferCoefficient",
        "convectiveHeatFlux",
        "liquidFrictionVelocity",
        "liquidTemperatureWallFunction",
        "quenchingHeatTransferCoefficient",
        "quenchingHeatFlux",
        "quenchingTemperatureProfile",
        "cellTemperatureProfile",
        "heatFluxBalanceTolerance",
        "maximumWallVaporFraction",
        "surfaceTension",
        "contactAngle",
        "bubbleNumberDensity",
        "coefficient",
        "vaporizationCoefficient",
        "nucleationVolumeFraction",
        "bubbleRadius",
        "mobility",
        "sourceLayers"
    };
    auto readCoefficientBlock = [&](const P* sourceBlock) {
        if (!sourceBlock) return;
        for (const std::string& key : coefficientKeys) {
            std::string value = jsonValueAfterKey(sourceBlock, key);
            if (!value.empty()) {
                config.coefficients[sourceTokenKey(key)] =
                    valueDoubleValue(value, 0.0);
            }
        }
    };
    readCoefficientBlock(block);
    readCoefficientBlock(modelSubBlock);

    const std::vector<std::string> selectionKeys = {
        "departureDiameterModel",
        "departureFrequencyModel",
        "nucleationSiteDensityModel",
        "criticalCavityRadiusModel",
        "liquidFrictionVelocityModel",
        "wallTemperatureModel",
        "convectiveHeatFluxModel",
        "quenchingHeatFluxModel",
        "heatFluxBudgetCheck",
        "sourcePatch",
        "sourcePatches",
        "wallPatch",
        "wallPatches",
        "patch",
        "patches",
        "sourceRegion",
        "liquidSpecies",
        "vaporSpecies"
    };
    auto readSelectionBlock = [&](const P* sourceBlock) {
        if (!sourceBlock) return;
        for (const std::string& key : selectionKeys) {
            std::string value = jsonValueAfterKey(sourceBlock, key);
            if (!value.empty()) {
                config.selections[sourceTokenKey(key)] = valueUnquote(value);
            }
        }
    };
    readSelectionBlock(block);
    readSelectionBlock(modelSubBlock);

    auto parseWallBoilingBlock = [&](const P* heatBlock,
                                     const std::string& fallbackPatch) {
        if (!heatBlock) return;
        std::string enabled = jsonValueAfterKey(heatBlock, "enabled");
        if (enabled.empty()) enabled = jsonValueAfterKey(heatBlock, "enable");
        if (!boolValue(enabled, true)) return;

        WallHeatSetting base;
        base.mode = "RPI";
        std::string heatFlux = jsonValueAfterKey(heatBlock, "qEvaporative");
        if (heatFlux.empty()) heatFlux = jsonValueAfterKey(heatBlock, "q");
        if (heatFlux.empty()) heatFlux = jsonValueAfterKey(heatBlock, "heatFlux");
        if (heatFlux.empty()) heatFlux = jsonValueAfterKey(heatBlock, "wallHeatFlux");
        if (heatFlux.empty()) heatFlux = jsonValueAfterKey(heatBlock, "value");
        if (!heatFlux.empty()) {
            base.heatFlux = valueDoubleValue(heatFlux, 0.0);
        } else {
            base.heatFlux = std::numeric_limits<double>::quiet_NaN();
        }
        std::string coordinate =
            valueLower(valueUnquote(
                jsonValueAfterKey(heatBlock, "rangeCoordinate")));
        if (coordinate.empty()) {
            coordinate = valueLower(valueUnquote(
                jsonValueAfterKey(heatBlock, "coordinate")));
        }
        if (!coordinate.empty()) {
            if (coordinate == "x") base.rangeCoordinate = 0;
            else if (coordinate == "y") base.rangeCoordinate = 1;
            else if (coordinate == "z") base.rangeCoordinate = 2;
            else fatalCaseConfig(
                "phaseChange wallBoiling coordinate supports x, y or z.");
            std::string range =
                jsonValueAfterKey(heatBlock, "range");
            if (range.empty()) {
                range = jsonValueAfterKey(heatBlock, "bounds");
            }
            const auto values = valueNumbers(range);
            if (values.size() != 2 || values[1] <= values[0]) {
                fatalCaseConfig(
                    "phaseChange wallBoiling range requires (min max).");
            }
            base.rangeMinimum = values[0];
            base.rangeMaximum = values[1];
        }
        std::string wallTemperature =
            jsonValueAfterKey(heatBlock, "wallTemperature");
        if (wallTemperature.empty()) {
            wallTemperature = jsonValueAfterKey(heatBlock, "Tw");
        }
        if (!wallTemperature.empty()) {
            base.wallTemperature =
                valueDoubleValue(wallTemperature, base.wallTemperature);
        }

        std::string patchText = jsonValueAfterKey(heatBlock, "patch");
        if (patchText.empty()) patchText = jsonValueAfterKey(heatBlock, "patches");
        if (patchText.empty()) patchText = jsonValueAfterKey(heatBlock, "set");
        if (patchText.empty()) patchText = jsonValueAfterKey(heatBlock, "zone");
        std::vector<std::string> patches = valueWordList(patchText);
        if (patches.empty() && !fallbackPatch.empty()) {
            patches.push_back(fallbackPatch);
        }
        if (patches.empty()) {
            fatalCaseConfig("phaseChange wallBoiling block must name a patch.");
        }
        for (const std::string& patch : patches) {
            WallHeatSetting setting = base;
            setting.patch = patch;
            setting.coupleEnergy = true;
            config.wallBoiling.push_back(setting);
        }
    };

    const P* wallBoiling = jsonBlockAfterKey(modelSubBlock, "wallBoiling");
    if (!wallBoiling) {
        wallBoiling = jsonBlockAfterKey(modelSubBlock, "wallBoilingPatches");
    }
    if (!wallBoiling) {
        wallBoiling = jsonBlockAfterKey(block, "wallBoiling");
    }
    if (!wallBoiling) {
        wallBoiling = jsonBlockAfterKey(block, "wallBoilingPatches");
    }
    if (wallBoiling) {
        const auto patchBlocks = jsonChildBlocks(wallBoiling);
        if (patchBlocks.empty()
            && (!jsonValueAfterKey(wallBoiling, "heatFlux").empty()
                || !jsonValueAfterKey(wallBoiling, "q").empty()
                || !jsonValueAfterKey(wallBoiling, "qEvaporative").empty())) {
            parseWallBoilingBlock(wallBoiling, "");
        }
        for (const auto& [name, heatBlock] : patchBlocks) {
            parseWallBoilingBlock(heatBlock, name);
        }
    } else if (!jsonValueAfterKey(block, "heatFlux").empty()
               || !jsonValueAfterKey(block, "q").empty()
               || !jsonValueAfterKey(block, "qEvaporative").empty()) {
        parseWallBoilingBlock(block, "");
    }

    multiPhaseConfig_.phaseChange = config;
    try {
        MP::validateMultiPhaseConfig(multiPhaseConfig_, context);
    } catch (const std::exception& e) {
        std::cerr << "[SF FATAL] " << e.what() << std::endl;
        std::exit(1);
    }
    multiPhaseLoaded_ = true;
}

void CaseAdapter::resolvePhaseChangeWallBoilingHeatFlux() {
    if (!multiPhaseLoaded_) return;
    auto& pc = multiPhaseConfig_.phaseChange;
    const std::string model = Physics::PhaseChange::normalizeModel(pc.model);
    if (!pc.enabled || (model != "rpi" && model != "wallboiling")) return;

    if (pc.wallBoiling.empty()) {
        for (const WallHeatSetting& heat : wallHeatBoundarySettings_) {
            if (!std::isfinite(heat.heatFlux) || heat.heatFlux <= 0.0) {
                fatalCaseConfig(
                    "T boundary wallHeatFlux patch '" + heat.patch
                    + "' used by RPI phaseChange must provide finite positive q.");
            }
            WallHeatSetting setting;
            setting.patch = heat.patch;
            setting.mode = "RPI";
            setting.heatFlux = heat.heatFlux;
            setting.coupleEnergy = true;
            pc.wallBoiling.push_back(setting);
        }
        if (pc.wallBoiling.empty()) {
            fatalCaseConfig(
                "RPI phaseChange requires at least one wallHeatFlux patch in "
                "the T boundary or an explicit wallBoiling patch list in "
                "models/phaseChange.yaml.");
        }
    }

    for (WallHeatSetting& setting : pc.wallBoiling) {
        if (std::isfinite(setting.heatFlux) && setting.heatFlux > 0.0) {
            continue;
        }

        const WallHeatSetting* match = nullptr;
        for (const WallHeatSetting& heat : wallHeatBoundarySettings_) {
            if (sourceTokenKey(heat.patch) != sourceTokenKey(setting.patch)) {
                continue;
            }
            if (!std::isfinite(heat.heatFlux) || heat.heatFlux <= 0.0) {
                fatalCaseConfig(
                    "T boundary wallHeatFlux patch '" + heat.patch
                    + "' used by RPI phaseChange must provide finite positive q.");
            }
            if (match != nullptr
                && std::abs(match->heatFlux - heat.heatFlux) > 1.0e-12) {
                fatalCaseConfig(
                    "RPI phaseChange patch '" + setting.patch
                    + "' matches multiple T boundary wallHeatFlux entries.");
            }
            match = &heat;
        }
        if (match == nullptr) {
            fatalCaseConfig(
                "RPI phaseChange wallBoiling patch '" + setting.patch
                + "' does not define heatFlux in models/phaseChange.yaml; "
                "provide a matching T boundary patch with type wallHeatFlux and q.");
        }
        setting.heatFlux = match->heatFlux;
    }

    try {
        Physics::Multiphase::validateMultiPhaseConfig(
            multiPhaseConfig_, "OpenFOAM phaseChange/T wallHeatFlux");
    } catch (const std::exception& e) {
        std::cerr << "[SF FATAL] " << e.what() << std::endl;
        std::exit(1);
    }
}

void CaseAdapter::decodeNativeSections() {
    decodeRuntimeSection();
    decodeOutputSection();
    decodeAlgorithmSection();
    if (!caseConfig_.createMesh) {
        // coupling preset 是 formulation/plan 输入；它不再携带 solver family。
        caseConfig_.solver.pressure = solverProperties_;
        caseConfig_.compatFlowLabel = legacyFlowLabel_;
    }

    // 模型的 presence 由 typed section 的 has* 决定，不再由磁盘文件名探测
    // 是否存在决定。顺序保持：物理/相模型先于数值块，数值块先于 field。
    decodeTurbulenceSection();
    decodePhaseSystemSection();
    decodePhaseChangeSection();
    decodeThermoDynamicsSection();
    decodeIBMSection();
    decodeILWSection();
    // Source models are decoded once, after registered field boundary contributions.
    // wallHeat boundary settings 由 field 解码阶段收集，必须延后处理。

    decodeNumericsSection();
    decodeParallelSection();

    if (!caseConfig_.createMesh) {
        decodeFields();
        resolvePhaseChangeWallBoilingHeatFlux();
    }

    scLoaded_ = true;
    if (!caseConfig_.createMesh) {
        icLoaded_ = true;
        bcLoaded_ = true;
    }
    FDM::BoundaryConfig& boundaries = caseConfig_.solver.boundaries;
    boundaries.ilwEnabled = boundaries.ilwEnabled && boundaries.ilwOrder > 0;
    if (!boundaries.ilwEnabled) boundaries.ilwOrder = 0;
    // ibmBoundary / ilwOrder 是 numerics 对边界闭合的读数；ILW 由 typed
    // boundaries 决定，而不是由字符串暂存值决定。
    caseConfig_.solver.numerics.ibmBoundary =
        boundaries.ilwEnabled ? FDM::IBMBoundaryScheme::ILW
                              : FDM::IBMBoundaryScheme::LowOrder;
    caseConfig_.solver.numerics.ilwOrder = boundaries.ilwOrder;
}

void CaseAdapter::decodePhaseSystemSection() {
    if (!sections_.hasPhaseSystem) return;
    const std::string context = "models/multiPhase.yaml";
    const P* text = &sections_.phaseSystem;

    namespace MP = Physics::Multiphase;
    MP::MultiPhaseConfig config;
    config.enabled = true;
    std::string value = jsonValueAfterKey(text, "multiPhases");
    if (value.empty()) value = jsonValueAfterKey(text, "multiphase");
    if (!value.empty()) config.enabled = valueBoolValue(value, config.enabled);

    const P* multiPhase = jsonBlockAfterKey(text, "multiPhase");
    if (!multiPhase) multiPhase = text;

    value = jsonValueAfterKey(multiPhase, "phaseSystem");
    if (value.empty()) value = jsonValueAfterKey(multiPhase, "type");
    if (value.empty()) value = jsonValueAfterKey(multiPhase, "model");
    if (!value.empty()) config.type = valueUnquote(value);

    if (!jsonValueAfterKey(multiPhase,"solver").empty()
        || !jsonValueAfterKey(multiPhase,"solverType").empty())
        throw std::runtime_error("Removed multiPhase solver-family selection: declare variables in STATE, equations in WHAT, execution in HOW and providers in WHICH.");

    value = jsonValueAfterKey(multiPhase, "defaultPhase");
    if (!value.empty()) config.defaultPhase = valueUnquote(value);
    value = jsonValueAfterKey(multiPhase, "defaultSignedDistance");
    if (value.empty()) value = jsonValueAfterKey(multiPhase, "signedDistance");
    if (!value.empty()) {
        config.defaultSignedDistance =
            valueDoubleValue(value, config.defaultSignedDistance);
    }
    value = jsonValueAfterKey(multiPhase, "initialPressure");
    if (value.empty()) value = jsonValueAfterKey(multiPhase, "p0");
    if (!value.empty()) config.initialPressure = valueDoubleValue(value, config.initialPressure);

    PhasePropertiesReader::parseFieldModels(
        text, multiPhase, config);

    PhasePropertiesReader::parsePhaseDefinitions(text, config);
    PhasePropertiesReader::parseEulerianProperties(
        text, multiPhase, config);
    PhasePropertiesReader::parseSetAssignments(text, config);

    try {
        MP::validateMultiPhaseConfig(config, context);
    } catch (const std::exception& e) {
        std::cerr << "[SF FATAL] " << e.what() << std::endl;
        std::exit(1);
    }

    multiPhaseConfig_ = config;
    multiPhaseLoaded_ = true;
}

std::optional<CaseConfig> CaseAdapter::decodeNativeCase() {
    broadcast("Loading case: ", caseFilePath_);
    caseConfig_ = {};
    // 每个语义块覆盖这份 typed 基线；解码过程完全不读 parser 全局暂存值。
    caseConfig_.solver = FDM::defaultSolverConfig();
    sourceScheme_.clear();

    // ── 1. native 语义块 → typed 规格 (不做字典翻译) ──
    decodeNativeSections();

    if (caseConfig_.createMesh) {
        // mesh generator 是 mesh.yaml 声明的格式适配器，不是 case 目录约定。
        if (!meshParameters_) {
            broadcast("Error: ",
                      "Missing mesh generator in mesh/mesh.yaml");
            return std::nullopt;
        }
    } else {
        for (const std::string& meshFile : meshFiles_) {
            const std::string meshPath = resolveCaseFile(caseDir_, meshFile);
            if (!fileExists(meshPath)) {
                broadcast("Error: ",
                          "Missing Sonic Fluid mesh file: " + meshPath);
                return std::nullopt;
            }
        }
    }

    if (!caseConfig_.createMesh) {
        decodeSourceSettings();
    }

    // VTK writer 保存独立快照，输出期间不再查询 parser 全局状态。
    densityOutputBC_ = caseConfig_.solver.boundaries.density;
    velocityOutputBC_ = caseConfig_.solver.boundaries.velocity;
    pressureOutputBC_ = caseConfig_.solver.boundaries.pressure;
    ibmOutputEnabled_ = caseConfig_.solver.ibm.enabled;

    // ── 3. 创建输出目录 ──
    if (outputDir_.empty()) outputDir_ = "result";
    std::string fullOut = caseDir_ + "/" + outputDir_;
    std::filesystem::create_directories(fullOut);
    broadcast("Output dir: ", fullOut);

    buildCaseConfig();

    return caseConfig_;
}

} // namespace SF
