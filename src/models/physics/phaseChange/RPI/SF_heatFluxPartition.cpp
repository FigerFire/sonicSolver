/// @file SF_heatFluxPartition.cpp
/// @brief RPI 壁面沸腾闭式模型与热流分配实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.07.07-----------*/

#include "SF_heatFluxPartition.h"
#include "SF_closure.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace SF {
namespace Physics {
namespace PhaseChange {
namespace RPI {

namespace {

std::string requireSelection(const Multiphase::PhaseChangeOptions& config,
                             const std::string& key,
                             const std::string& context) {
    const std::string raw = selection(config, key, "");
    if (raw.empty()) {
        throw std::runtime_error(
            context + ": RPI phaseChange requires explicit " + key
            + " selection.");
    }
    return normalizeChoice(raw);
}

double requireNonNegative(const Multiphase::PhaseChangeOptions& config,
                          const std::string& key,
                          const std::string& context) {
    const double value = coefficient(config, key, -1.0);
    if (!std::isfinite(value) || value < 0.0) {
        throw std::runtime_error(
            context + ": RPI phaseChange requires " + key
            + " to be finite and >= 0.");
    }
    return value;
}

void validateSensibleModel(const Multiphase::PhaseChangeOptions& config,
                           const std::string& modelKey,
                           const std::string& hKey,
                           const std::string& qKey,
                           const std::string& context) {
    const std::string mode = requireSelection(config, modelKey, context);
    if (mode == "none" || mode == "off" || mode == "disabled") return;
    if (mode == "papertable1" || mode == "table1") {
        if (modelKey == "convectiveHeatFluxModel") {
            const std::string frictionModel = normalizeChoice(
                selection(
                    config, "liquidFrictionVelocityModel", "explicit"));
            if (frictionModel == "explicit"
                || frictionModel == "constant") {
                (void)requireNonNegative(
                    config, "liquidFrictionVelocity", context);
            } else if (frictionModel != "localwallshear") {
                throw std::runtime_error(
                    context + ": liquidFrictionVelocityModel supports "
                    "explicit or localWallShear.");
            }
            (void)requireNonNegative(
                config, "liquidTemperatureWallFunction", context);
        } else {
            (void)requireNonNegative(config, hKey, context);
            (void)requireNonNegative(
                config, "quenchingTemperatureProfile", context);
            (void)requireNonNegative(
                config, "cellTemperatureProfile", context);
        }
        return;
    }
    if (mode == "coefficient" || mode == "heattransfercoefficient") {
        (void)requireNonNegative(config, hKey, context);
        return;
    }
    if (mode == "fixed" || mode == "explicit") {
        (void)requireNonNegative(config, qKey, context);
        return;
    }
    throw std::runtime_error(
        context + ": unsupported RPI " + modelKey + " '"
        + selection(config, modelKey, "") + "'. Supported values are "
        "none, paperTable1, coefficient and fixed.");
}

double sensibleHeatFlux(const ModelContext& ctx,
                        const std::string& modelKey,
                        const std::string& hKey,
                        const std::string& qKey,
                        double wallTemperature,
                        double liquidTemperature) {
    const auto& config = ctx.config.phaseChange;
    const std::string mode = normalizeChoice(selection(config, modelKey, ""));
    if (mode == "none" || mode == "off" || mode == "disabled") return 0.0;
    if (mode == "fixed" || mode == "explicit") {
        return requireNonNegative(config, qKey, "RPI phaseChange");
    }
    const double deltaT = wallTemperature - liquidTemperature;
    if (!std::isfinite(deltaT)) {
        throw std::runtime_error(
            "RPI phaseChange received non-finite wall-liquid temperature difference.");
    }
    double q = 0.0;
    if (mode == "papertable1" || mode == "table1") {
        if (modelKey == "convectiveHeatFluxModel") {
            const std::string frictionModel = normalizeChoice(
                selection(
                    config, "liquidFrictionVelocityModel", "explicit"));
            const double frictionVelocity =
                frictionModel == "localwallshear"
                ? ctx.localFrictionVelocity
                : requireNonNegative(
                    config, "liquidFrictionVelocity", "RPI phaseChange");
            const double temperatureWallFunction = requireNonNegative(
                config, "liquidTemperatureWallFunction", "RPI phaseChange");
            if (!std::isfinite(frictionVelocity)
                || frictionVelocity < 0.0
                || temperatureWallFunction <= 0.0
                || !std::isfinite(ctx.liquidPhase.density)
                || ctx.liquidPhase.density <= 0.0
                || !std::isfinite(ctx.liquidPhase.specificHeat)
                || ctx.liquidPhase.specificHeat <= 0.0) {
                throw std::runtime_error(
                    "RPI Table-1 convection needs positive rho_l, Cp_l "
                    "and liquidTemperatureWallFunction.");
            }
            q = ctx.liquidPhase.density * ctx.liquidPhase.specificHeat
                * frictionVelocity / temperatureWallFunction * deltaT;
        } else {
            const double profileAtQuench = requireNonNegative(
                config, "quenchingTemperatureProfile", "RPI phaseChange");
            const double profileAtCell = requireNonNegative(
                config, "cellTemperatureProfile", "RPI phaseChange");
            if (profileAtCell <= 0.0) {
                throw std::runtime_error(
                    "RPI Table-1 quenching needs positive "
                    "cellTemperatureProfile.");
            }
            q = requireNonNegative(config, hKey, "RPI phaseChange")
                * profileAtQuench / profileAtCell * deltaT;
        }
    } else {
        q = requireNonNegative(config, hKey, "RPI phaseChange") * deltaT;
    }
    if (!std::isfinite(q) || q < 0.0) {
        throw std::runtime_error(
            "RPI phaseChange sensible heat partition became negative or non-finite.");
    }
    return q;
}

} // namespace

void validateHeatFluxConfig(const Multiphase::PhaseChangeOptions& config,
                            const std::string& context) {
    validateSensibleModel(config, "convectiveHeatFluxModel",
                          "convectiveHeatTransferCoefficient",
                          "convectiveHeatFlux", context);
    validateSensibleModel(config, "quenchingHeatFluxModel",
                          "quenchingHeatTransferCoefficient",
                          "quenchingHeatFlux", context);
    const std::string wallTemperatureModel = normalizeChoice(
        selection(config, "wallTemperatureModel", "explicit"));
    if (wallTemperatureModel != "explicit"
        && wallTemperatureModel != "fixed"
        && wallTemperatureModel != "heatfluxbalance") {
        throw std::runtime_error(
            context + ": wallTemperatureModel supports explicit or "
            "heatFluxBalance.");
    }
}

bool budgetCheckEnabled(const Multiphase::PhaseChangeOptions& config) {
    const std::string raw = normalizeChoice(
        selection(config, "heatFluxBudgetCheck", "true"));
    return !(raw == "false" || raw == "off" || raw == "no" || raw == "0");
}

HeatFluxPartition evaluateHeatFluxPartition(const ModelContext& ctx,
                                            double wallTemperature,
                                            double liquidTemperature,
                                            double wallMassFlux) {
    if (!std::isfinite(wallTemperature) || wallTemperature <= 0.0
        || !std::isfinite(liquidTemperature) || liquidTemperature <= 0.0) {
        throw std::runtime_error(
            "RPI phaseChange invalid wall/near-wall temperature.");
    }
    if (!std::isfinite(wallMassFlux) || wallMassFlux <= 0.0) {
        throw std::runtime_error(
            "RPI phaseChange wall mass flux must be finite and > 0.");
    }

    const auto& pc = ctx.config.phaseChange;
    HeatFluxPartition q;
    q.convective = sensibleHeatFlux(ctx, "convectiveHeatFluxModel",
                                    "convectiveHeatTransferCoefficient",
                                    "convectiveHeatFlux",
                                    wallTemperature,
                                    liquidTemperature);
    q.quenching = sensibleHeatFlux(ctx, "quenchingHeatFluxModel",
                                   "quenchingHeatTransferCoefficient",
                                   "quenchingHeatFlux",
                                   wallTemperature,
                                   liquidTemperature);
    q.evaporative = wallMassFlux * pc.latentHeat;
    q.total = q.convective + q.quenching + q.evaporative;
    if (!std::isfinite(q.total) || q.total < 0.0) {
        throw std::runtime_error(
            "RPI phaseChange heat-flux partition is non-finite or negative.");
    }
    return q;
}

WallHeatBalance evaluateWallHeatBalance(
        const ModelContext& ctx,
        double liquidTemperature,
        double referenceWallTemperature,
        double wallHeatFlux) {
    if (!std::isfinite(wallHeatFlux) || wallHeatFlux <= 0.0) {
        throw std::runtime_error(
            "RPI wall heat balance requires finite positive wallHeatFlux.");
    }
    const auto evaluate = [&](double wallTemperature) {
        WallHeatBalance state;
        state.wallTemperature = wallTemperature;
        state.closure = evaluateClosure(
            ctx, wallTemperature, liquidTemperature);
        state.heat = evaluateHeatFluxPartition(
            ctx, wallTemperature, liquidTemperature,
            state.closure.wallMassFlux);
        return state;
    };

    const std::string model = normalizeChoice(selection(
        ctx.config.phaseChange, "wallTemperatureModel", "explicit"));
    if (model == "explicit" || model == "fixed") {
        if (!std::isfinite(referenceWallTemperature)
            || referenceWallTemperature <= 0.0) {
            throw std::runtime_error(
                "RPI explicit wallTemperatureModel requires wallTemperature.");
        }
        return evaluate(referenceWallTemperature);
    }
    if (model != "heatfluxbalance") {
        throw std::runtime_error(
            "RPI wallTemperatureModel supports explicit or heatFluxBalance.");
    }

    const double saturationTemperature =
        ctx.config.phaseChange.saturationTemperature;
    double lower = std::nextafter(
        std::max(saturationTemperature, liquidTemperature),
        std::numeric_limits<double>::max());
    WallHeatBalance lowerState = evaluate(lower);
    if (lowerState.heat.total >= wallHeatFlux) {
        throw std::runtime_error(
            "RPI heatFluxBalance has no nucleate-boiling root: the minimum "
            "partition already exceeds qWall.");
    }

    double upper = referenceWallTemperature;
    if (!std::isfinite(upper) || upper <= lower) {
        upper = lower + 1.0;
    }
    WallHeatBalance upperState = evaluate(upper);
    double span = upper - lower;
    for (int expansion = 0;
         upperState.heat.total < wallHeatFlux && expansion < 32;
         ++expansion) {
        span *= 2.0;
        upper = lower + span;
        upperState = evaluate(upper);
    }
    if (upperState.heat.total < wallHeatFlux) {
        throw std::runtime_error(
            "RPI heatFluxBalance could not bracket the wall-temperature root.");
    }

    WallHeatBalance midpointState = upperState;
    for (int iteration = 0; iteration < 80; ++iteration) {
        const double midpoint = 0.5 * (lower + upper);
        midpointState = evaluate(midpoint);
        const double residual = midpointState.heat.total - wallHeatFlux;
        if (std::abs(residual) <= 1.0e-10 * wallHeatFlux) {
            return midpointState;
        }
        if (residual < 0.0) {
            lower = midpoint;
            lowerState = midpointState;
        } else {
            upper = midpoint;
            upperState = midpointState;
        }
    }
    throw std::runtime_error(
        "RPI heatFluxBalance wall-temperature solve did not converge.");
}

} // namespace RPI
} // namespace PhaseChange
} // namespace Physics
} // namespace SF
