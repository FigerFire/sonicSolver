/// @file SF_closure.cpp
/// @brief RPI 壁面沸腾闭式模型与热流分配实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.07.07-----------*/

#include "SF_closure.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace SF {
namespace Physics {
namespace PhaseChange {
namespace RPI {

namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

double requirePositive(const Multiphase::PhaseChangeOptions& config,
                       const std::string& key,
                       const std::string& context) {
    const double value = coefficient(config, key, 0.0);
    if (!std::isfinite(value) || value <= 0.0) {
        throw std::runtime_error(
            context + ": RPI phaseChange requires " + key
            + " to be finite and > 0.");
    }
    return value;
}

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

const Multiphase::PhaseProperties& vaporPhase(const ModelContext& ctx) {
    return ctx.alphaIsLiquid ? ctx.otherPhase : ctx.alphaPhase;
}

const Multiphase::PhaseProperties& liquidPhase(const ModelContext& ctx) {
    return ctx.liquidPhase;
}

double departureDiameter(
        const ModelContext& ctx,
        double wallTemperature,
        double liquidTemperature) {
    const auto& pc = ctx.config.phaseChange;
    const std::string model =
        requireSelection(pc, "departureDiameterModel", "RPI phaseChange");
    if (model == "explicit" || model == "constant") {
        return requirePositive(
            pc, "departureDiameter", "RPI phaseChange");
    }
    const double subcooling =
        pc.saturationTemperature - liquidTemperature;
    if (!std::isfinite(subcooling) || subcooling <= 0.0) {
        throw std::runtime_error(
            "RPI departure-diameter correlation requires positive liquid "
            "subcooling.");
    }
    if (model == "tolubinskykostanchuk" || model == "tolubinsky") {
        const double diameter0 =
            requirePositive(pc, "referenceDepartureDiameter", "RPI");
        const double subcooling0 =
            requirePositive(pc, "referenceSubcooling", "RPI");
        return diameter0 * std::exp(-subcooling / subcooling0);
    }
    if (model == "unal") {
        const double pressure = std::isfinite(ctx.localPressure)
            ? ctx.localPressure : pc.saturationPressure;
        const double a = requirePositive(pc, "unalA", "RPI");
        const double velocity = std::isfinite(ctx.localLiquidSpeed)
            ? ctx.localLiquidSpeed
            : requirePositive(pc, "nearWallLiquidVelocity", "RPI");
        const double velocity0 =
            requirePositive(pc, "unalReferenceVelocity", "RPI");
        const auto& liquid = liquidPhase(ctx);
        const auto& vapor = vaporPhase(ctx);
        const double b = subcooling
            / (2.0 * (1.0 - vapor.density / liquid.density));
        const double phi = velocity >= velocity0
            ? std::pow(velocity / velocity0, 0.47) : 1.0;
        if (!std::isfinite(pressure) || pressure <= 0.0
            || !std::isfinite(b) || b <= 0.0
            || !std::isfinite(phi) || phi <= 0.0) {
            throw std::runtime_error(
                "RPI Unal correlation received invalid p, b or phi.");
        }
        (void)wallTemperature;
        return 2.42e-5 * std::pow(pressure, 0.709) * a
             / std::sqrt(b * phi);
    }
    throw std::runtime_error(
        "RPI departureDiameterModel supports explicit, "
        "TolubinskyKostanchuk or Unal.");
}

double departureFrequency(
        const ModelContext& ctx,
        double diameter) {
    const auto& pc = ctx.config.phaseChange;
    const std::string model =
        requireSelection(pc, "departureFrequencyModel", "RPI phaseChange");
    if (model == "explicit" || model == "constant") {
        return requirePositive(
            pc, "departureFrequency", "RPI phaseChange");
    }
    if (model == "cole") {
        const double gravity =
            requirePositive(pc, "gravityMagnitude", "RPI");
        const auto& liquid = liquidPhase(ctx);
        const auto& vapor = vaporPhase(ctx);
        const double argument = 4.0 * gravity
            * (liquid.density - vapor.density)
            / (3.0 * diameter * liquid.density);
        if (!std::isfinite(argument) || argument <= 0.0) {
            throw std::runtime_error(
                "RPI Cole frequency produced an invalid square-root argument.");
        }
        return std::sqrt(argument);
    }
    throw std::runtime_error(
        "RPI departureFrequencyModel supports explicit or Cole.");
}

double nucleationDensity(
        const ModelContext& ctx,
        double wallTemperature,
        double liquidTemperature) {
    const auto& pc = ctx.config.phaseChange;
    const std::string model = requireSelection(
        pc, "nucleationSiteDensityModel", "RPI phaseChange");
    if (model == "explicit" || model == "constant") {
        return requirePositive(
            pc, "nucleationSiteDensity", "RPI phaseChange");
    }
    const double superheat =
        wallTemperature - pc.saturationTemperature;
    if (!std::isfinite(superheat) || superheat <= 0.0) {
        throw std::runtime_error(
            "RPI nucleation correlation requires positive wall superheat.");
    }
    if (model == "lemmertchawla") {
        const double calibration =
            requirePositive(pc, "lemmertChawlaM", "RPI");
        const double exponent =
            requirePositive(pc, "lemmertChawlaExponent", "RPI");
        return std::pow(calibration * superheat, exponent);
    }
    if (model == "hibikiishii") {
        const double averageDensity =
            requirePositive(pc, "averageCavityDensity", "RPI");
        const double contactAngle =
            requirePositive(pc, "contactAngle", "RPI");
        const double angleScale =
            requirePositive(pc, "contactAngleScale", "RPI");
        const double cavityLength =
            requirePositive(pc, "cavityLengthScale", "RPI");
        const double densityFunction =
            requirePositive(pc, "hibikiDensityFunction", "RPI");
        double criticalRadius = 0.0;
        const std::string radiusModel = normalizeChoice(
            selection(pc, "criticalCavityRadiusModel", "explicit"));
        if (radiusModel == "explicit" || radiusModel == "constant") {
            criticalRadius =
                requirePositive(pc, "criticalCavityRadius", "RPI");
        } else if (radiusModel == "localthermo") {
            const double surfaceTension =
                Multiphase::phasePairSurfaceTension(
                    ctx.config, liquidPhase(ctx).name,
                    vaporPhase(ctx).name);
            const double superheat =
                wallTemperature - pc.saturationTemperature;
            const auto& vapor = vaporPhase(ctx);
            criticalRadius =
                2.0 * surfaceTension * pc.saturationTemperature
                / (vapor.density * pc.latentHeat * superheat);
        } else {
            throw std::runtime_error(
                "RPI criticalCavityRadiusModel supports explicit or "
                "localThermo.");
        }
        (void)liquidTemperature;
        const double wetting = 1.0 - std::exp(
            -contactAngle * contactAngle
            / (8.0 * angleScale * angleScale));
        const double cavity = std::exp(
            densityFunction * cavityLength / criticalRadius) - 1.0;
        const double result = averageDensity * wetting * cavity;
        if (!std::isfinite(result) || result <= 0.0) {
            throw std::runtime_error(
                "RPI HibikiIshii correlation produced invalid site density.");
        }
        return result;
    }
    throw std::runtime_error(
        "RPI nucleationSiteDensityModel supports explicit, "
        "LemmertChawla or HibikiIshii.");
}

} // namespace

std::string normalizeChoice(std::string value) {
    std::string out;
    out.reserve(value.size());
    for (char c : value) {
        if (c == '_' || c == '-' || std::isspace((unsigned char)c)) continue;
        out.push_back((char)std::tolower((unsigned char)c));
    }
    return out;
}

void validateClosureConfig(const Multiphase::PhaseChangeOptions& config,
                           const std::string& context) {
    const std::string diameter =
        requireSelection(config, "departureDiameterModel", context);
    const std::string frequency =
        requireSelection(config, "departureFrequencyModel", context);
    const std::string nucleation =
        requireSelection(config, "nucleationSiteDensityModel", context);
    if (diameter == "explicit" || diameter == "constant") {
        (void)requirePositive(config, "departureDiameter", context);
    } else if (diameter == "tolubinskykostanchuk"
               || diameter == "tolubinsky") {
        (void)requirePositive(config, "referenceDepartureDiameter", context);
        (void)requirePositive(config, "referenceSubcooling", context);
    } else if (diameter == "unal") {
        (void)requirePositive(config, "unalA", context);
        (void)requirePositive(config, "nearWallLiquidVelocity", context);
        (void)requirePositive(config, "unalReferenceVelocity", context);
    } else {
        throw std::runtime_error(
            context + ": unsupported RPI departureDiameterModel.");
    }
    if (frequency == "explicit" || frequency == "constant") {
        (void)requirePositive(config, "departureFrequency", context);
    } else if (frequency == "cole") {
        (void)requirePositive(config, "gravityMagnitude", context);
    } else {
        throw std::runtime_error(
            context + ": unsupported RPI departureFrequencyModel.");
    }
    if (nucleation == "explicit" || nucleation == "constant") {
        (void)requirePositive(config, "nucleationSiteDensity", context);
    } else if (nucleation == "lemmertchawla") {
        (void)requirePositive(config, "lemmertChawlaM", context);
        (void)requirePositive(config, "lemmertChawlaExponent", context);
    } else if (nucleation == "hibikiishii") {
        for (const std::string& key : {
                 "averageCavityDensity", "contactAngle",
                 "contactAngleScale", "cavityLengthScale",
                 "hibikiDensityFunction"}) {
            (void)requirePositive(config, key, context);
        }
        const std::string radiusModel = normalizeChoice(
            selection(config, "criticalCavityRadiusModel", "explicit"));
        if (radiusModel == "explicit" || radiusModel == "constant") {
            (void)requirePositive(config, "criticalCavityRadius", context);
        } else if (radiusModel == "localthermo") {
            (void)requirePositive(config, "surfaceTension", context);
        } else {
            throw std::runtime_error(
                context + ": unsupported RPI criticalCavityRadiusModel.");
        }
    } else {
        throw std::runtime_error(
            context + ": unsupported RPI nucleationSiteDensityModel.");
    }
}

BubbleClosure evaluateClosure(const ModelContext& ctx,
                              double wallTemperature,
                              double liquidTemperature) {
    if (!std::isfinite(wallTemperature) || wallTemperature <= 0.0) {
        throw std::runtime_error(
            "RPI phaseChange received invalid wall temperature.");
    }
    const auto& pc = ctx.config.phaseChange;
    const auto& vapor = vaporPhase(ctx);
    if (!std::isfinite(vapor.density) || vapor.density <= 0.0) {
        throw std::runtime_error(
            "RPI phaseChange requires vapor phase density > 0.");
    }

    BubbleClosure closure;
    closure.departureDiameter = departureDiameter(
        ctx, wallTemperature, liquidTemperature);
    closure.departureFrequency =
        departureFrequency(ctx, closure.departureDiameter);
    closure.nucleationSiteDensity =
        nucleationDensity(ctx, wallTemperature, liquidTemperature);

    const double bubbleVolume =
        pi * closure.departureDiameter * closure.departureDiameter
        * closure.departureDiameter / 6.0;
    closure.wallMassFlux =
        closure.nucleationSiteDensity * bubbleVolume * vapor.density
        * closure.departureFrequency;
    if (!std::isfinite(closure.wallMassFlux)
        || closure.wallMassFlux <= 0.0) {
        throw std::runtime_error(
            "RPI phaseChange produced invalid wall mass flux from d_d, f_d "
            "and N_a.");
    }
    return closure;
}

} // namespace RPI
} // namespace PhaseChange
} // namespace Physics
} // namespace SF
