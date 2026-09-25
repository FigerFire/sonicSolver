/// @file SF_phaseChange.cpp
/// @brief 相变模型、饱和性质与守恒传递实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.07.07-----------*/

#include "SF_phaseChange.h"

#include "RPI/SF_RPI.h"
#include "SF_scalarTransport.h"
#include "SF_utility.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace SF {
namespace Physics {
namespace PhaseChange {

namespace {

std::string compactLower(std::string value) {
    std::string out;
    out.reserve(value.size());
    for (char c : value) {
        if (c == '_' || c == '-' || std::isspace((unsigned char)c)) continue;
        out.push_back((char)std::tolower((unsigned char)c));
    }
    return out;
}

const Multiphase::PhaseProperties* findByRole(
        const Multiphase::MultiPhaseConfig& config,
        Multiphase::PhaseRole role) {
    for (const auto& phase : config.phases) {
        const Multiphase::PhaseRole resolved =
            phase.role == Multiphase::PhaseRole::Unknown
                ? Multiphase::parsePhaseRole(phase.name) : phase.role;
        if (resolved == role) return &phase;
    }
    return nullptr;
}

const Multiphase::PhaseProperties* findAlphaPhase(
        const Multiphase::MultiPhaseConfig& config) {
    if (const auto* phase =
            Multiphase::findPhase(config, config.alpha.phaseName)) {
        return phase;
    }
    const Multiphase::PhaseRole role =
        Multiphase::parsePhaseRole(config.alpha.phaseName);
    if (role != Multiphase::PhaseRole::Unknown) {
        return findByRole(config, role);
    }
    return nullptr;
}

Multiphase::PhaseRole resolvedRole(
        const Multiphase::PhaseProperties& phase) {
    return phase.role == Multiphase::PhaseRole::Unknown
        ? Multiphase::parsePhaseRole(phase.name) : phase.role;
}

const Multiphase::PhaseProperties* findOtherPhase(
        const Multiphase::MultiPhaseConfig& config,
        const Multiphase::PhaseProperties* alphaPhase) {
    if (alphaPhase == nullptr) return nullptr;
    const std::string alphaName =
        Multiphase::normalizePhaseName(alphaPhase->name);
    for (const auto& phase : config.phases) {
        if (Multiphase::normalizePhaseName(phase.name) != alphaName) {
            return &phase;
        }
    }
    return nullptr;
}

const Multiphase::PhaseProperties* findLiquidPhase(
        const Multiphase::MultiPhaseConfig& config,
        const Multiphase::PhaseProperties* alphaPhase,
        const Multiphase::PhaseProperties* otherPhase) {
    if (alphaPhase && resolvedRole(*alphaPhase)
        == Multiphase::PhaseRole::Liquid) {
        return alphaPhase;
    }
    if (otherPhase && resolvedRole(*otherPhase)
        == Multiphase::PhaseRole::Liquid) {
        return otherPhase;
    }
    return findByRole(config, Multiphase::PhaseRole::Liquid);
}

void requirePhase(const Multiphase::PhaseProperties* phase,
                  const std::string& label) {
    if (phase == nullptr) {
        throw std::runtime_error("PhaseChange: missing " + label + " phase.");
    }
    if (!std::isfinite(phase->density) || phase->density <= 0.0) {
        throw std::runtime_error(
            "PhaseChange: phase '" + phase->name
            + "' density must be finite and > 0.");
    }
    if (!std::isfinite(phase->specificHeat) || phase->specificHeat <= 0.0) {
        throw std::runtime_error(
            "PhaseChange: phase '" + phase->name
            + "' specificHeat must be finite and > 0.");
    }
}

void dispatchRates(const ModelContext& ctx, std::vector<double>& mdot) {
    const std::string model = normalizeModel(ctx.config.phaseChange.model);
    if (model == "lee" || model == "superheat") {
        addLeeRates(ctx, mdot);
    } else if (model == "stefan") {
        addStefanRates(ctx, mdot);
    } else if (model == "hks" || model == "schrage"
               || model == "hertzknudsenschrage") {
        addHKSRates(ctx, mdot);
    } else if (model == "enthalpyporosity" || model == "mushyzone") {
        addEnthalpyPorosityRates(ctx, mdot);
    } else if (model == "rpi" || model == "wallboiling") {
        RPI::computeRates(ctx, mdot);
    } else if (model == "schnerrsauer") {
        addSchnerrSauerRates(ctx, mdot);
    } else if (model == "zgb" || model == "zwartgerberbelamri") {
        addZGBRates(ctx, mdot);
    } else if (model == "phasefield" || model == "allen-cahn"
               || model == "allencahn") {
        addPhaseFieldRates(ctx, mdot);
    } else if (model == "saturationproperties"
               || model == "thermophysicalclosure") {
        addSaturationPropertyRates(ctx, mdot);
    } else {
        throw std::runtime_error(
            "PhaseChange: unsupported model '" + ctx.config.phaseChange.model
            + "'. Supported models: Lee, Stefan, HKS, EnthalpyPorosity, RPI, "
              "SchnerrSauer, ZGB, PhaseField, SaturationProperties.");
    }
}

/// @brief 输出聚合诊断日志。
void logDiagnostics(const PhaseChangeDiagnostics& diag) {
    if (diag.totalFluidCells == 0) return;

    std::ostringstream oss;
    oss << "active=" << diag.activeCells
        << " of " << diag.totalFluidCells;
    if (diag.activeCells == 0) {
        SF::broadcast("PhaseChange diag  : ", oss.str());
        return;
    }
    if (diag.limitedCells > 0) {
        oss << ", limited=" << diag.limitedCells
            << " (mass=" << diag.failedMassCells
            << " heat=" << diag.failedHeatCells << ")";
    }
    oss << ", alpha=[" << diag.minAlpha << ", " << diag.maxAlpha << "]";
    oss << ", T=[" << diag.minTemperature << ", " << diag.maxTemperature
        << "]";
    oss << ", mdot=[" << diag.minMdot << ", " << diag.maxMdot << "]";
    oss << ", S_E=[" << diag.minEnergySource << ", "
        << diag.maxEnergySource << "]";
    if (diag.limitedCells > 0) {
        oss << ", maxRelReduction=" << diag.maxRateReduction;
    }
    if (diag.failedPostStateCells > 0) {
        oss << ", FAILED_POST_STATE=" << diag.failedPostStateCells;
    }
    if (diag.firstBadCellIndex >= 0) {
        oss << ", firstBadCell=" << diag.firstBadCellIndex;
    }
    SF::broadcast("PhaseChange diag  : ", oss.str());
}

} // namespace

// ---- 公共入口 ----

std::string normalizeModel(std::string model) {
    std::string m = compactLower(std::move(model));
    if (m == "pri") return "rpi";
    return m;
}

double coefficient(const Multiphase::PhaseChangeOptions& config,
                   const std::string& key,
                   double fallback) {
    auto it = config.coefficients.find(compactLower(key));
    return it == config.coefficients.end() ? fallback : it->second;
}

std::string selection(const Multiphase::PhaseChangeOptions& config,
                      const std::string& key,
                      const std::string& fallback) {
    auto it = config.selections.find(compactLower(key));
    return it == config.selections.end() ? fallback : it->second;
}

void validateConfig(const Multiphase::MultiPhaseConfig& config,
                    const std::string& context) {
    const auto& pc = config.phaseChange;
    const std::string model = normalizeModel(pc.model);
    if (!pc.enabled || model == "none") return;
    if (!pc.phaseNames.empty()) {
        if (pc.phaseNames.size() != 2
            || Multiphase::normalizePhaseName(pc.phaseNames[0])
                == Multiphase::normalizePhaseName(pc.phaseNames[1])
            || Multiphase::findPhase(config, pc.phaseNames[0]) == nullptr
            || Multiphase::findPhase(config, pc.phaseNames[1]) == nullptr) {
            throw std::runtime_error(
                context + ": phaseChange phases must name two distinct "
                "declared phases.");
        }
    } else if (config.phases.size() != 2) {
        throw std::runtime_error(
            context + ": phaseChange in a system with more than two phases "
            "must declare phases (source target).");
    }

    const bool perPhaseTemperature =
        Multiphase::isEulerianEulerianType(config.type)
        && std::all_of(config.phases.begin(), config.phases.end(),
            [](const Multiphase::PhaseProperties& phase) {
                return phase.initialTemperatureDeclared;
            });
    if (!config.temperature.enabled && !perPhaseTemperature) {
        throw std::runtime_error(
            context + ": [phaseChange] requires mixture temperature or "
            "Eulerian per-phase temperature fields.");
    }
    if (!pc.energyCoupling) {
        throw std::runtime_error(
            context + ": [phaseChange].energyCoupling cannot be false; "
            "phase change must consume or release latent heat.");
    }
    if (!config.alpha.bounded
        || !std::isfinite(config.alpha.lowerBound)
        || !std::isfinite(config.alpha.upperBound)
        || config.alpha.lowerBound > config.alpha.upperBound) {
        throw std::runtime_error(
            context + ": [phaseChange] requires finite ordered [alpha] bounds "
            "so mass-transfer rates can be conservatively limited.");
    }
    if (!std::isfinite(pc.latentHeat) || pc.latentHeat <= 0.0) {
        throw std::runtime_error(
            context + ": [phaseChange].latentHeat must be finite and > 0.");
    }
    for (const auto& phase : config.phases) {
        if (!std::isfinite(phase.specificHeat) || phase.specificHeat <= 0.0) {
            throw std::runtime_error(
                context + ": phase '" + phase.name
                + "' needs finite positive specificHeat when phaseChange is enabled.");
        }
    }

    const bool thermalModel =
        model == "lee" || model == "superheat" || model == "stefan"
        || model == "hks" || model == "enthalpyporosity"
        || model == "rpi" || model == "phasefield";
    if (thermalModel
        && (!std::isfinite(pc.saturationTemperature)
            || pc.saturationTemperature <= 0.0)) {
        throw std::runtime_error(
            context + ": [phaseChange].saturationTemperature/Tsat must be "
            "finite and > 0 for model '" + pc.model + "'.");
    }
    if ((model == "schnerrsauer" || model == "zgb")
        && (!std::isfinite(pc.saturationPressure)
            || pc.saturationPressure <= 0.0)) {
        throw std::runtime_error(
            context + ": [phaseChange].saturationPressure/psat must be "
            "finite and > 0 for cavitation models.");
    }
    if ((model == "lee" || model == "superheat")
        && pc.evaporationCoefficient == 0.0
        && pc.condensationCoefficient == 0.0) {
        throw std::runtime_error(
            context + ": Lee/superheat phaseChange needs a non-zero "
            "evaporation or condensation coefficient.");
    }
    if (model == "rpi" || model == "wallboiling") {
        RPI::validateConfig(pc, context);
    }
    if (model == "saturationproperties"
        || model == "thermophysicalclosure") {
        throw std::runtime_error(
            context + ": SaturationProperties is a property closure, not an "
            "independent mass-transfer source. Choose Lee, Stefan, HKS, RPI, "
            "SchnerrSauer, ZGB, EnthalpyPorosity, or PhaseField.");
    }
}

PhaseChangeRateResult computeRates(const Field& field,
                                   const ScalarField& alpha,
                                   const ScalarField& temperature,
                                   const Multiphase::MultiPhaseConfig& config,
                                   double dt) {
    if (!alpha.isCompatibleWith(field) || !temperature.isCompatibleWith(field)) {
        throw std::runtime_error(
            "PhaseChange: alpha/temperature dimensions do not match Field.");
    }
    if (!std::isfinite(dt) || dt <= 0.0) {
        throw std::runtime_error("PhaseChange: dt must be finite and > 0.");
    }

    const auto* alphaPhase = findAlphaPhase(config);
    const auto* otherPhase = findOtherPhase(config, alphaPhase);
    if (!config.phaseChange.phaseNames.empty()) {
        const auto* first = Multiphase::findPhase(
            config, config.phaseChange.phaseNames[0]);
        const auto* second = Multiphase::findPhase(
            config, config.phaseChange.phaseNames[1]);
        const std::string alphaName = alphaPhase
            ? Multiphase::normalizePhaseName(alphaPhase->name) : "";
        if (first
            && Multiphase::normalizePhaseName(first->name) == alphaName) {
            otherPhase = second;
        } else if (second
                   && Multiphase::normalizePhaseName(second->name)
                       == alphaName) {
            otherPhase = first;
        } else {
            throw std::runtime_error(
                "PhaseChange: transported alpha phase is not part of the "
                "configured phaseChange phases.");
        }
    }
    const auto* liquidPhase =
        findLiquidPhase(config, alphaPhase, otherPhase);
    requirePhase(alphaPhase, "alpha");
    requirePhase(otherPhase, "secondary");
    requirePhase(liquidPhase, "liquid");
    const auto alphaRole = resolvedRole(*alphaPhase);
    if (alphaRole != Multiphase::PhaseRole::Liquid
        && alphaRole != Multiphase::PhaseRole::Gas) {
        throw std::runtime_error(
            "PhaseChange: alpha phase must have liquid or gas role.");
    }

    const bool alphaIsLiquid = alphaRole == Multiphase::PhaseRole::Liquid;

    PhaseChangeRateResult result;
    result.mdot.assign((size_t)field.TotalSize(), 0.0);

    // Rate law 读取前先验证输入；任何模型内部的端点正则化都不得掩盖
    // 非有限温度或越界 alpha。
    Math::forInterior(field, [&](int i, int j, int k) {
        if (field.CellFlag(i, j, k) != FLUID_CELL) return;
        const double a = alpha(i,j,k);
        const double T = temperature(i,j,k);
        if (!std::isfinite(a) || a < 0.0 || a > 1.0
            || !std::isfinite(T) || T <= 0.0) {
            std::ostringstream oss;
            oss << "PhaseChange: invalid rate-law input at ("
                << i << "," << j << "," << k << "): alpha=" << a
                << ", T=" << T << ".";
            throw std::runtime_error(oss.str());
        }
    });

    ModelContext ctx{field, alpha, temperature, config, *alphaPhase,
                     *otherPhase, *liquidPhase, alphaIsLiquid, dt};
    dispatchRates(ctx, result.mdot);

    Math::forInterior(field, [&](int i, int j, int k) {
        if (field.CellFlag(i, j, k) != FLUID_CELL) return;
        ++result.diagnostics.totalFluidCells;

        const int id = alpha.getIdx(i, j, k);
        const double rate = result.mdot[(size_t)id];
        if (!std::isfinite(rate)) {
            std::ostringstream oss;
            oss << "PhaseChange: non-finite mass source at ("
                << i << "," << j << "," << k << ").";
            throw std::runtime_error(oss.str());
        }
        if (rate == 0.0) return;

        ++result.diagnostics.activeCells;
        const double oldAlpha = alpha(i, j, k);
        const double oldT = temperature(i, j, k);
        if (!std::isfinite(oldAlpha) || !std::isfinite(oldT) || oldT <= 0.0) {
            std::ostringstream oss;
            oss << "PhaseChange: invalid alpha/T at ("
                << i << "," << j << "," << k << ").";
            throw std::runtime_error(oss.str());
        }

        if (oldAlpha < 0.0 || oldAlpha > 1.0) {
            std::ostringstream oss;
            oss << "PhaseChange: alpha outside [0,1] at ("
                << i << "," << j << "," << k << "): " << oldAlpha
                << ". Rate laws do not clamp transported state.";
            throw std::runtime_error(oss.str());
        }

        // ---- 聚合统计 ----
        result.diagnostics.minAlpha = std::min(result.diagnostics.minAlpha, oldAlpha);
        result.diagnostics.maxAlpha = std::max(result.diagnostics.maxAlpha, oldAlpha);
        result.diagnostics.minTemperature = std::min(result.diagnostics.minTemperature, oldT);
        result.diagnostics.maxTemperature = std::max(result.diagnostics.maxTemperature, oldT);
        result.diagnostics.minMdot = std::min(result.diagnostics.minMdot, rate);
        result.diagnostics.maxMdot = std::max(result.diagnostics.maxMdot, rate);
        result.diagnostics.minEnergySource = 0.0;
        result.diagnostics.maxEnergySource = 0.0;
    });

    return result;
}

void reportDiagnostics(const PhaseChangeDiagnostics& diagnostics) {
    logDiagnostics(diagnostics);
}

} // namespace PhaseChange
} // namespace Physics
} // namespace SF
