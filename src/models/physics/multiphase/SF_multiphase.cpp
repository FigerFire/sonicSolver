/// @file SF_multiphase.cpp
/// @brief 多相配置、物性与模型一致性校验实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_multiphase.h"

#include "SF_scalarTransport.h"
#include "SF_utility.h"
#include "SF_phaseChange.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace SF {
namespace Physics {
namespace Multiphase {

namespace {
const PhaseProperties* findByRole(const MultiPhaseConfig& config,
                                  PhaseRole role) {
    for (const PhaseProperties& phase : config.phases) {
        const PhaseRole resolved =
            phase.role == PhaseRole::Unknown
                ? parsePhaseRole(phase.name) : phase.role;
        if (resolved == role) return &phase;
    }
    return nullptr;
}

const PhaseProperties* findAlphaPhase(const MultiPhaseConfig& config) {
    if (const PhaseProperties* phase = findPhase(config, config.alpha.phaseName)) {
        return phase;
    }
    const PhaseRole role = parsePhaseRole(config.alpha.phaseName);
    if (role != PhaseRole::Unknown) return findByRole(config, role);
    return nullptr;
}

const PhaseProperties* findOtherPhase(const MultiPhaseConfig& config,
                                      const PhaseProperties* alphaPhase) {
    if (alphaPhase == nullptr) return nullptr;

    const PhaseRole alphaRole =
        alphaPhase->role == PhaseRole::Unknown
            ? parsePhaseRole(alphaPhase->name) : alphaPhase->role;
    if (alphaRole == PhaseRole::Liquid) {
        if (const PhaseProperties* gas = findByRole(config, PhaseRole::Gas)) {
            return gas;
        }
    }
    if (alphaRole == PhaseRole::Gas) {
        if (const PhaseProperties* liquid = findByRole(config, PhaseRole::Liquid)) {
            return liquid;
        }
    }

    const std::string alphaName = normalizePhaseName(alphaPhase->name);
    for (const PhaseProperties& phase : config.phases) {
        if (normalizePhaseName(phase.name) != alphaName) return &phase;
    }
    return nullptr;
}

PhaseRole resolvedPhaseRole(const PhaseProperties& phase) {
    return phase.role == PhaseRole::Unknown
        ? parsePhaseRole(phase.name) : phase.role;
}

const PhaseProperties* findLiquidPhase(const MultiPhaseConfig& config,
                                       const PhaseProperties* alphaPhase,
                                       const PhaseProperties* otherPhase) {
    if (alphaPhase && resolvedPhaseRole(*alphaPhase) == PhaseRole::Liquid) {
        return alphaPhase;
    }
    if (otherPhase && resolvedPhaseRole(*otherPhase) == PhaseRole::Liquid) {
        return otherPhase;
    }
    return findByRole(config, PhaseRole::Liquid);
}

void requireMixturePhase(const PhaseProperties* phase,
                         const std::string& label) {
    if (phase == nullptr) {
        throw std::runtime_error(
            "Mixture alpha properties: missing " + label + " phase.");
    }
    if (!std::isfinite(phase->density) || phase->density <= 0.0) {
        throw std::runtime_error(
            "Mixture alpha properties: phase '" + phase->name
            + "' density must be finite and > 0.");
    }
    if (!std::isfinite(phase->viscosity) || phase->viscosity < 0.0) {
        throw std::runtime_error(
            "Mixture alpha properties: phase '" + phase->name
            + "' viscosity must be finite and >= 0.");
    }
}

FDM::Scalar::Bounds alphaBounds(const MultiPhaseConfig& config) {
    FDM::Scalar::Bounds bounds;
    bounds.enabled = config.alpha.bounded;
    bounds.lower = config.alpha.lowerBound;
    bounds.upper = config.alpha.upperBound;
    bounds.mode = FDM::Scalar::parseBoundMode(config.alpha.boundMode);
    return bounds;
}

FDM::Scalar::TransportConfig alphaTransportConfig(
        const MultiPhaseConfig& config) {
    FDM::Scalar::TransportConfig transport;
    transport.name = config.alpha.fieldName;
    transport.convectionEnabled = config.alpha.transportEnabled;
    transport.diffusionEnabled = config.alpha.diffusionEnabled;
    transport.diffusivity = config.alpha.diffusivity;
    transport.bounds = alphaBounds(config);
    return transport;
}

double trackedMixturePhaseDensity(const MultiPhaseConfig& config) {
    const PhaseProperties* phase = findAlphaPhase(config);
    requireMixturePhase(phase, "alpha");
    return phase->density;
}

std::vector<BCSetting<double>> phaseMassBoundaries(
        const MultiPhaseConfig& config) {
    std::vector<BCSetting<double>> result = config.alpha.boundaryConditions;
    const double rhoPhase = trackedMixturePhaseDensity(config);
    for (BCSetting<double>& bc : result) {
        if (bc.type == FIXED_VALUE) bc.value *= rhoPhase;
    }
    return result;
}

FDM::Scalar::TransportConfig phaseMassTransportConfig(
        const MultiPhaseConfig& config) {
    FDM::Scalar::TransportConfig transport = alphaTransportConfig(config);
    const double rhoPhase = trackedMixturePhaseDensity(config);
    transport.name = "phaseMass." + config.alpha.phaseName;
    if (transport.bounds.enabled) {
        transport.bounds.lower *= rhoPhase;
        transport.bounds.upper *= rhoPhase;
    }
    return transport;
}

FDM::Scalar::Bounds temperatureBounds(const MultiPhaseConfig& config) {
    FDM::Scalar::Bounds bounds;
    bounds.enabled = config.temperature.bounded;
    bounds.lower = config.temperature.lowerBound;
    bounds.upper = config.temperature.upperBound;
    bounds.mode = FDM::Scalar::parseBoundMode(config.temperature.boundMode);
    return bounds;
}

FDM::Scalar::TransportConfig temperatureTransportConfig(
        const MultiPhaseConfig& config) {
    FDM::Scalar::TransportConfig transport;
    transport.name = config.temperature.fieldName;
    transport.convectionEnabled = config.temperature.transportEnabled;
    transport.diffusionEnabled = config.temperature.diffusionEnabled;
    transport.diffusivity = config.temperature.diffusivity;
    transport.bounds = temperatureBounds(config);
    return transport;
}

bool phaseChangeActive(const MultiPhaseConfig& config) {
    return config.phaseChange.enabled
        && normalizeModelType(config.phaseChange.model) != "none";
}

/// @brief 返回显式定容比热；缺失时拒绝热力学闭合。
/// @param phase 相属性。
/// @return Cv 或 specificHeat。
double effectiveCv(const PhaseProperties& phase) {
    if (!std::isfinite(phase.Cv) || phase.Cv <= 0.0) {
        throw std::runtime_error(
            "Thermodynamic closure: phase '" + phase.name
            + "' requires explicit finite Cv > 0; Cp is not an internal-energy fallback.");
    }
    return phase.Cv;
}

/// @brief 返回有效参考内能：e0 若已显式设置，否则为 0。
/// @param phase 相属性。
/// @return e0 或 0.0。
double effectiveE0(const PhaseProperties& phase) {
    return std::isfinite(phase.e0) ? phase.e0 : 0.0;
}

double mixtureRhoCv(double alpha,
                    const PhaseProperties& alphaPhase,
                    const PhaseProperties& otherPhase) {
    const double cvAlpha = effectiveCv(alphaPhase);
    const double cvOther = effectiveCv(otherPhase);
    const double rhoCv =
        alpha * alphaPhase.density * cvAlpha
        + (1.0 - alpha) * otherPhase.density * cvOther;
    if (!std::isfinite(rhoCv) || rhoCv <= 0.0) {
        throw std::runtime_error(
            "Mixture temperature source: rho*Cv must be finite and > 0."
            " (alpha=" + std::to_string(alpha)
            + " cvAlpha=" + std::to_string(cvAlpha)
            + " cvOther=" + std::to_string(cvOther) + ")");
    }
    return rhoCv;
}

double projectedAlphaValue(double value,
                           const FDM::Scalar::Bounds& bounds) {
    if (!bounds.enabled || bounds.mode != FDM::Scalar::BoundMode::Project) {
        return value;
    }
    return std::max(bounds.lower, std::min(bounds.upper, value));
}

double cellVolume(const Field& field, int i, int j, int k) {
    const double inverseJacobian = field.Jac(i, j, k);
    if (!std::isfinite(inverseJacobian)
        || std::abs(inverseJacobian) <= 1.0e-300) {
        std::ostringstream oss;
        oss << "MultiPhase diagnostics: invalid inverse Jacobian at ("
            << i << "," << j << "," << k << ") value="
            << inverseJacobian;
        throw std::runtime_error(oss.str());
    }
    return std::abs(1.0 / inverseJacobian);
}

double relativeError(double value, double reference) {
    const double scale = std::max(std::abs(reference), 1.0e-300);
    return value / scale;
}

double flowMassIntegral(const Field& field) {
    double sum = 0.0;
    Math::forInterior(field, [&](int i, int j, int k) {
        if (field.CellFlag(i, j, k) != FLUID_CELL) return;
        const double rho = field(i, j, k, RHO);
        if (!std::isfinite(rho) || rho <= 0.0) {
            std::ostringstream oss;
            oss << "MultiPhase diagnostics: invalid flow density at ("
                << i << "," << j << "," << k << ") value=" << rho;
            throw std::runtime_error(oss.str());
        }
        sum += rho * cellVolume(field, i, j, k);
    });
    return sum;
}

double trackedPhaseDensity(const MultiPhaseConfig& config) {
    const PhaseProperties* phase = findAlphaPhase(config);
    requireMixturePhase(phase, "alpha");
    return phase->density;
}

} // namespace

MultiPhaseModel::MultiPhaseModel() = default;
MultiPhaseModel::~MultiPhaseModel() = default;
MultiPhaseModel::MultiPhaseModel(MultiPhaseModel&&) noexcept = default;
MultiPhaseModel& MultiPhaseModel::operator=(MultiPhaseModel&&) noexcept = default;

const ConservationDiagnostics& MultiPhaseModel::conservation() const {
    return diagnostics_;
}

void MultiPhaseModel::configure(const MultiPhaseConfig& config) {
    validateMultiPhaseConfig(config, "multiPhase");
    config_ = config;
    initialized_ = false;
    hasMixtureProperties_ = false;
    hasTemperature_ = false;
    hasPhaseChangeSource_ = false;
    diagnostics_ = ConservationDiagnostics();
    if (isLevelSetType(config_.type)) {
        throw std::runtime_error(
            "MultiPhaseModel no longer owns Level Set; create an "
            "InterfaceModels::Model through makeInterfaceModel().");
    }
}

void MultiPhaseModel::initialize(const Field& field) {
    if (isThermal()) {
        initializeTemperature(field);
        initialized_ = true;
        return;
    }
    if (isMixture()) {
        initializeAlpha(field);
        return;
    }
    throw std::runtime_error(
        "MultiPhaseModel: unsupported model type '" + config_.type + "'.");
}

void MultiPhaseModel::initializeAlpha(const Field& field) {
    if (!isMixture()) {
        throw std::runtime_error(
            "MultiPhaseModel: initializeAlpha requires type = mixture.");
    }
    if (!std::isfinite(config_.alpha.defaultValue)) {
        throw std::runtime_error(
            "MultiPhaseModel: alpha defaultValue must be finite.");
    }
    alpha_.setupLike(field, config_.alpha.fieldName,
                     config_.alpha.defaultValue);
    FDM::Scalar::applyInitialConditions(
        field, alpha_, config_.alpha.initialConditions);
    FDM::Scalar::applyBoundaryConditions(
        field, alpha_, config_.alpha.boundaryConditions,
        scalarILWEnabled_, scalarILWAccuracyOrder_);
    FDM::Scalar::applyBounds(field, alpha_, alphaBounds(config_));
    phaseMass_.setupLike(field,
                         "phaseMass." + config_.alpha.phaseName,
                         0.0);
    updatePhaseMassFromAlpha(field);
    if (config_.temperature.enabled) {
        initializeTemperature(field);
    } else {
        hasTemperature_ = false;
    }
    phaseChangeEnergySource_.assign((size_t)field.TotalSize(), 0.0);
    hasPhaseChangeSource_ = false;
    updateMixtureProperties(field);
    initialized_ = true;
    updateConservationDiagnostics(field, true);
}

void MultiPhaseModel::advance(Field& field, double dt) {
    if (!enabled()) return;
    if (!initialized_) {
        throw std::runtime_error(
            "MultiPhaseModel::advance called before model initialization.");
    }
    if (!std::isfinite(dt) || dt <= 0.0) {
        throw std::runtime_error(
            "MultiPhaseModel::advance requires a finite positive dt.");
    }

    if (isThermal()) {
        advanceTemperature(field, dt);
        return;
    }

    if (isMixture()) {
        throw std::runtime_error(
            "MultiPhaseModel::advance cannot update mixture state outside "
            "the solver integrator; use assembleIntegratedRHS().");
    }

    throw std::runtime_error(
        "MultiPhaseModel::advance supports only legacy thermal state; "
        "mixture is advanced by assembleIntegratedRHS().");
}

void MultiPhaseModel::updateMixtureProperties(const Field& field) {
    if (!alpha_.isCompatibleWith(field)) {
        throw std::runtime_error(
            "Mixture alpha properties: alpha dimensions do not match Field.");
    }

    const PhaseProperties* alphaPhase = findAlphaPhase(config_);
    const PhaseProperties* otherPhase =
        findOtherPhase(config_, alphaPhase);
    requireMixturePhase(alphaPhase, "alpha");
    requireMixturePhase(otherPhase, "secondary");

    mixtureDensity_.assign((size_t)field.TotalSize(), 0.0);
    mixtureViscosity_.assign((size_t)field.TotalSize(), 0.0);

    for (int k = 0; k < field.MZ(); ++k) {
        for (int j = 0; j < field.MY(); ++j) {
            for (int i = 0; i < field.MX(); ++i) {
                const int id = alpha_.getIdx(i, j, k);
                const double a = alpha_(i, j, k);
                if (!std::isfinite(a)) {
                    std::ostringstream oss;
                    oss << "Mixture alpha properties: non-finite alpha at ("
                        << i << "," << j << "," << k << ").";
                    throw std::runtime_error(oss.str());
                }
                mixtureDensity_[(size_t)id] =
                    a * alphaPhase->density
                    + (1.0 - a) * otherPhase->density;
                mixtureViscosity_[(size_t)id] =
                    a * alphaPhase->viscosity
                    + (1.0 - a) * otherPhase->viscosity;
            }
        }
    }

    hasMixtureProperties_ = true;
}

void MultiPhaseModel::initializeTemperature(const Field& field) {
    temperature_.setupLike(field,
                           config_.temperature.fieldName,
                           config_.temperature.defaultValue);
    FDM::Scalar::applyInitialConditions(
        field, temperature_, config_.temperature.initialConditions);
    FDM::Scalar::applyBoundaryConditions(
        field, temperature_, config_.temperature.boundaryConditions,
        scalarILWEnabled_, scalarILWAccuracyOrder_);
    FDM::Scalar::applyBounds(field, temperature_, temperatureBounds(config_));
    Math::forInterior(field, [&](int i, int j, int k) {
        if (field.CellFlag(i, j, k) != FLUID_CELL) return;
        const double T = temperature_(i, j, k);
        if (!std::isfinite(T) || T <= 0.0) {
            std::ostringstream oss;
            oss << "Mixture temperature: cell (" << i << "," << j << ","
                << k << ") is not covered by an explicit positive "
                << "temperature initial condition.";
            throw std::runtime_error(oss.str());
        }
    });
    hasTemperature_ = true;
}

void MultiPhaseModel::advanceTemperature(Field& field, double dt) {
    if (!config_.temperature.enabled) return;
    if (!hasTemperature_) initializeTemperature(field);

    if (config_.temperature.transportEnabled
        || config_.temperature.diffusionEnabled) {
        auto transport = temperatureTransportConfig(config_);
        transport.ilwBoundaryEnabled = scalarILWEnabled_;
        transport.ilwAccuracyOrder = scalarILWAccuracyOrder_;
        FDM::Scalar::advanceForwardEuler(
            field,
            temperature_,
            config_.temperature.boundaryConditions,
            transport,
            dt);
    } else {
        FDM::Scalar::applyBoundaryConditions(
            field, temperature_, config_.temperature.boundaryConditions,
            scalarILWEnabled_, scalarILWAccuracyOrder_);
        FDM::Scalar::applyBounds(
            field, temperature_, temperatureBounds(config_));
    }
}


void MultiPhaseModel::refreshDerivedFromConserved(const Field& field) {
    if (!hasTemperature_) {
        throw std::runtime_error(
            "deriveTemperature: temperature field is not initialized.");
    }
    if (!hasMixtureProperties_) {
        updateMixtureProperties(field);
    }

    const PhaseProperties* alphaPhase = findAlphaPhase(config_);
    const PhaseProperties* otherPhase = findOtherPhase(config_, alphaPhase);
    requireMixturePhase(alphaPhase, "alpha");
    requireMixturePhase(otherPhase, "secondary");

    const double cvAlpha    = effectiveCv(*alphaPhase);
    const double cvOther    = effectiveCv(*otherPhase);
    const double e0Alpha    = effectiveE0(*alphaPhase);
    const double e0Other    = effectiveE0(*otherPhase);
    const double rhoAlpha   = alphaPhase->density;
    const double rhoOther   = otherPhase->density;
    const double T_ref_alpha = std::isfinite(alphaPhase->T_ref)
        ? alphaPhase->T_ref : 0.0;
    const double T_ref_other = std::isfinite(otherPhase->T_ref)
        ? otherPhase->T_ref : 0.0;

    Math::forInterior(field, [&](int i, int j, int k) {
        if (field.CellFlag(i, j, k) != FLUID_CELL) return;

        const double rho  = field(i, j, k, RHO);
        const double ru   = field(i, j, k, RU);
        const double rv   = field(i, j, k, RV);
        const double rw   = field(i, j, k, RW);
        const double rhoE = field(i, j, k, E);

        if (!std::isfinite(rho) || rho <= 0.0) {
            throw std::runtime_error(
                "deriveTemperature: non-finite or non-positive rho.");
        }

        const double ke = 0.5 * (ru * ru + rv * rv + rw * rw) / rho;
        const double rhoEint = rhoE - ke;

        const double alpha   = alpha_(i, j, k);
        const double alpha_v = 1.0 - alpha;

        const double rhoCv = alpha * rhoAlpha * cvAlpha
                              + alpha_v * rhoOther * cvOther;
        if (!std::isfinite(rhoCv) || rhoCv <= 0.0) {
            throw std::runtime_error(
                "deriveTemperature: non-finite or non-positive rho*Cv.");
        }

        const double rhoEnergyOffset =
            alpha * rhoAlpha * (e0Alpha - cvAlpha * T_ref_alpha)
            + alpha_v * rhoOther * (e0Other - cvOther * T_ref_other);
        const double T = (rhoEint - rhoEnergyOffset) / rhoCv;
        if (!std::isfinite(T) || T <= 0.0) {
            std::ostringstream oss;
            throw std::runtime_error(oss.str());
        }

        temperature_(i, j, k) = T;
    });
}

void MultiPhaseModel::updatePhaseMassFromAlpha(const Field& field) {
    const double rhoPhase = trackedMixturePhaseDensity(config_);
    for (int k = 0; k < field.MZ(); ++k)
        for (int j = 0; j < field.MY(); ++j)
            for (int i = 0; i < field.MX(); ++i) {
                const double a = alpha_(i, j, k);
                if (!std::isfinite(a)) {
                    throw std::runtime_error(
                        "phase-mass initialization found non-finite alpha.");
                }
                phaseMass_(i, j, k) = rhoPhase * a;
            }
}

void MultiPhaseModel::updateAlphaFromPhaseMass(const Field& field) {
    const double rhoPhase = trackedMixturePhaseDensity(config_);
    for (int k = 0; k < field.MZ(); ++k)
        for (int j = 0; j < field.MY(); ++j)
            for (int i = 0; i < field.MX(); ++i) {
                const double mass = phaseMass_(i, j, k);
                if (!std::isfinite(mass)) {
                    throw std::runtime_error(
                        "phase-mass closure found non-finite conserved mass.");
                }
                alpha_(i, j, k) = mass / rhoPhase;
            }
}

void MultiPhaseModel::initializeConservedFromPrimitive(Field& field) {
    if (!isMixture() || !hasTemperature_) return;
    updateAlphaFromPhaseMass(field);
    const PhaseProperties* alphaPhase = findAlphaPhase(config_);
    const PhaseProperties* otherPhase = findOtherPhase(config_, alphaPhase);
    requireMixturePhase(alphaPhase, "alpha");
    requireMixturePhase(otherPhase, "secondary");
    const double cvAlpha = effectiveCv(*alphaPhase);
    const double cvOther = effectiveCv(*otherPhase);
    const double e0Alpha = effectiveE0(*alphaPhase);
    const double e0Other = effectiveE0(*otherPhase);
    Math::forInterior(field, [&](int i, int j, int k) {
        if (field.CellFlag(i, j, k) != FLUID_CELL) return;
        const double a = alpha_(i, j, k);
        const double T = temperature_(i, j, k);
        const double rho = a * alphaPhase->density
                         + (1.0 - a) * otherPhase->density;
        const double oldRho = field(i, j, k, RHO);
        if (!std::isfinite(rho) || rho <= 0.0
            || !std::isfinite(oldRho) || oldRho <= 0.0
            || !std::isfinite(T) || T <= 0.0) {
            throw std::runtime_error(
                "cold-start thermo closure found invalid alpha/T/rho.");
        }
        const double u = field(i, j, k, RU) / oldRho;
        const double v = field(i, j, k, RV) / oldRho;
        const double w = field(i, j, k, RW) / oldRho;
        const double inputPressure = Numerics::pressure(field, i, j, k);
        field(i, j, k, RHO) = rho;
        field(i, j, k, RU) = rho * u;
        field(i, j, k, RV) = rho * v;
        field(i, j, k, RW) = rho * w;
        const double internal =
            a * alphaPhase->density
                * (e0Alpha + cvAlpha * (T - alphaPhase->T_ref))
            + (1.0 - a) * otherPhase->density
                * (e0Other + cvOther * (T - otherPhase->T_ref));
        field(i, j, k, E) = internal
            + 0.5 * rho * (u*u + v*v + w*w);
        const double closedPressure = Numerics::pressure(field, i, j, k);
        const double pressureScale = std::max(std::abs(inputPressure), 1.0);
        if (!std::isfinite(inputPressure) || inputPressure <= 0.0
            || !std::isfinite(closedPressure) || closedPressure <= 0.0
            || std::abs(closedPressure - inputPressure)
                > 1.0e-8 * pressureScale) {
            std::ostringstream oss;
            oss << "cold-start thermo inputs are inconsistent at ("
                << i << "," << j << "," << k << "): input p="
                << inputPressure << ", EOS p(alpha,T)=" << closedPressure
                << ". Configure explicit Cv/e0/T_ref for every phase; "
                   "rhoE cannot satisfy independent p and T definitions.";
            throw std::runtime_error(oss.str());
        }
    });
    refreshDerivedFromConserved(field);
}

void MultiPhaseModel::assembleIntegratedRHS(
        Field& field, Residual& residual,
        std::vector<double>& phaseMassRHS,
        double dt) {
    if (!isMixture()) return;
    applyAuxiliaryBoundaryConditions(field);
    FDM::Scalar::assembleRHS(
        field, phaseMass_, phaseMassTransportConfig(config_), phaseMassRHS);
    if (!phaseChangeActive(config_)) {
        return;
    }
    if (!hasTemperature_) {
        throw std::runtime_error(
            "assembleIntegratedRHS: temperature field is not initialized.");
    }

    // 1. Thermo closure: derive T from conserved rhoE + alpha
    refreshDerivedFromConserved(field);

    // 2. Compute fresh phase change sources from current stage state
    const PhaseChange::PhaseChangeRateResult result =
        PhaseChange::computeRates(
            field, alpha_, temperature_, config_, dt);
    phaseChangeDiagnostics_ = result.diagnostics;

    // 3. Write S_E directly to Field::Source(E)
    //    Also cache for VTK/diagnostics (no-op for source term)
    phaseChangeEnergySource_.assign((size_t)field.TotalSize(), 0.0);
    Math::forInterior(field, [&](int i, int j, int k) {
        if (field.CellFlag(i, j, k) != FLUID_CELL) return;
        const int id = alpha_.getIdx(i, j, k);
        constexpr double SE = 0.0;
        phaseChangeEnergySource_[(size_t)id] = SE;
        residual.source(i, j, k, E) += SE;
    });
    hasPhaseChangeSource_ = false; // already written to Source, not needed later

    // 4. 被跟踪相的守恒质量 RHS；总能量不添加潜热源。
    const auto* alphaPhase = findAlphaPhase(config_);
    const auto* otherPhase = findOtherPhase(config_, alphaPhase);
    const auto* liquidPhase = otherPhase;
    // find the liquid phase
    if (alphaPhase && resolvedPhaseRole(*alphaPhase) == PhaseRole::Liquid)
        liquidPhase = alphaPhase;
    
    if (!alphaPhase || !liquidPhase) {
        throw std::runtime_error(
            "phase-change RHS cannot resolve tracked and liquid phases.");
    }
    const auto alphaRole = resolvedPhaseRole(*alphaPhase);
    const bool alphaIsLiquid = alphaRole == PhaseRole::Liquid;
    Math::forInterior(field, [&](int i, int j, int k) {
        if (field.CellFlag(i, j, k) != FLUID_CELL) return;
        const int id = alpha_.getIdx(i, j, k);
        const double mdot = result.mdot[(size_t)id];
        phaseMassRHS[(size_t)id] += alphaIsLiquid ? -mdot : mdot;
    });
}

void MultiPhaseModel::commitTimeLevel(const Field& field) {
    if (!enabled() || !initialized()) return;
    if (isMixture()) {
        applyAuxiliaryBoundaryConditions(field);
        refreshDerivedFromConserved(field);
        updateMixtureProperties(field);
        updateConservationDiagnostics(field, false);
        if (phaseChangeActive(config_)) {
            PhaseChange::reportDiagnostics(phaseChangeDiagnostics_);
        }
        return;
    }
    updateFlowCoupling(field);
}

void MultiPhaseModel::applyAuxiliaryBoundaryConditions(const Field& field) {
    if (!enabled()) return;
    if (!initialized()) {
        throw std::runtime_error(
            "MultiPhaseModel::applyAuxiliaryBoundaryConditions called before initialization.");
    }

    if (isThermal()) {
        if (config_.temperature.enabled && hasTemperature_) {
            FDM::Scalar::applyBoundaryConditions(
                field, temperature_, config_.temperature.boundaryConditions,
                scalarILWEnabled_, scalarILWAccuracyOrder_);
            FDM::Scalar::applyBounds(
                field, temperature_, temperatureBounds(config_));
        }
        return;
    }

    if (isMixture()) {
        const auto boundaries = phaseMassBoundaries(config_);
        const auto transport = phaseMassTransportConfig(config_);
        FDM::Scalar::applyBoundaryConditions(
            field, phaseMass_, boundaries,
            scalarILWEnabled_, scalarILWAccuracyOrder_);
        FDM::Scalar::applyBounds(field, phaseMass_, transport.bounds);
        updateAlphaFromPhaseMass(field);
        return;
    }

    throw std::runtime_error(
        "MultiPhaseModel::applyAuxiliaryBoundaryConditions: unsupported model type '"
        + config_.type + "'.");
}

void MultiPhaseModel::updateFlowCoupling(const Field& field) {
    if (!enabled()) return;
    if (!initialized()) {
        throw std::runtime_error(
            "MultiPhaseModel::updateFlowCoupling called before initialization.");
    }

    applyAuxiliaryBoundaryConditions(field);
    if (isThermal()) {
        return;
    }

    if (isMixture()) {
        updateAlphaFromPhaseMass(field);
        updateMixtureProperties(field);
        return;
    }

    throw std::runtime_error(
        "MultiPhaseModel::updateFlowCoupling: unsupported legacy model type '"
        + config_.type + "'.");
}

void MultiPhaseModel::addSourceTerms(Field& field, Residual& residual) const {
    if (!enabled()) return;
    if (!initialized()) {
        throw std::runtime_error(
            "MultiPhaseModel::addSourceTerms called before initialization.");
    }

    if (isThermal()) {
        return;
    }

    if (isMixture() && hasPhaseChangeSource_) {
        if ((int)phaseChangeEnergySource_.size() != field.TotalSize()) {
            throw std::runtime_error(
                "MultiPhase phaseChange source: source cache size does not match Field.");
        }
        Math::forInterior(field, [&](int i, int j, int k) {
            if (field.CellFlag(i, j, k) != FLUID_CELL) return;
            const int id = alpha_.getIdx(i, j, k);
            const double source = phaseChangeEnergySource_[(size_t)id];
            if (!std::isfinite(source)) {
                std::ostringstream oss;
                oss << "MultiPhase phaseChange source: non-finite energy "
                    << "source at (" << i << "," << j << "," << k << ").";
                throw std::runtime_error(oss.str());
            }
            residual.source(i, j, k, E) += source;
        });
    }
}

void MultiPhaseModel::applyBoundary(const Field& field) {
    applyAuxiliaryBoundaryConditions(field);
}

void MultiPhaseModel::correct(const Field& field, double dt) {
    (void)dt;
    updateFlowCoupling(field);
}

double MultiPhaseModel::dynamicViscosity(const Field& field,
                                         int i,
                                         int j,
                                         int k,
                                         double laminarMu) const {
    if (!enabled()) return laminarMu;
    if (!initialized()) {
        throw std::runtime_error(
            "MultiPhaseModel::dynamicViscosity called before initialization.");
    }

    double value = laminarMu;
    if (isMixture()) {
        if (!hasMixtureProperties_
            || (int)mixtureViscosity_.size() != field.TotalSize()) {
            throw std::runtime_error(
                "MultiPhaseModel::dynamicViscosity: mixture viscosity cache is not current.");
        }
        value = mixtureViscosity_[(size_t)alpha_.getIdx(i, j, k)];
    } else if (isThermal()) {
        value = laminarMu;
    } else {
        throw std::runtime_error(
            "MultiPhaseModel::dynamicViscosity: unsupported model type '"
            + config_.type + "'.");
    }

    if (!std::isfinite(value) || value < 0.0) {
        std::ostringstream oss;
        oss << "MultiPhaseModel::dynamicViscosity: invalid viscosity at ("
            << i << "," << j << "," << k << ") value=" << value;
        throw std::runtime_error(oss.str());
    }
    return value;
}

void MultiPhaseModel::updateConservationDiagnostics(
        const Field& field,
        bool resetReference) {
    if (!enabled()) return;

    if (!isMixture()) {
        throw std::runtime_error(
            "MultiPhaseModel diagnostics now own only legacy mixture state.");
    }
    const double phaseDensity = trackedPhaseDensity(config_);
    double phaseIndicator = 0.0;
    Math::forInterior(field, [&](int i, int j, int k) {
        if (field.CellFlag(i, j, k) != FLUID_CELL) return;
        const double indicator = alpha_(i, j, k);
        if (!std::isfinite(indicator)) {
            std::ostringstream oss;
            oss << "MultiPhase diagnostics: non-finite phase indicator at ("
                << i << "," << j << "," << k << ").";
            throw std::runtime_error(oss.str());
        }
        phaseIndicator += indicator * cellVolume(field, i, j, k);
    });

    const double phaseMass = phaseIndicator * phaseDensity;
    const double flowMass = flowMassIntegral(field);

    if (resetReference || !diagnostics_.initialized) {
        diagnostics_.initialized = true;
        diagnostics_.referencePhaseIndicator = phaseIndicator;
        diagnostics_.referencePhaseMass = phaseMass;
        diagnostics_.referenceFlowMass = flowMass;
    }

    diagnostics_.currentPhaseIndicator = phaseIndicator;
    diagnostics_.phaseIndicatorError =
        phaseIndicator - diagnostics_.referencePhaseIndicator;
    diagnostics_.phaseIndicatorRelativeError =
        relativeError(diagnostics_.phaseIndicatorError,
                      diagnostics_.referencePhaseIndicator);
    diagnostics_.currentPhaseMass = phaseMass;
    diagnostics_.phaseMassError =
        phaseMass - diagnostics_.referencePhaseMass;
    diagnostics_.phaseMassRelativeError =
        relativeError(diagnostics_.phaseMassError,
                      diagnostics_.referencePhaseMass);
    diagnostics_.currentFlowMass = flowMass;
    diagnostics_.flowMassError =
        flowMass - diagnostics_.referenceFlowMass;
    diagnostics_.flowMassRelativeError =
        relativeError(diagnostics_.flowMassError,
                      diagnostics_.referenceFlowMass);
}

} // namespace Multiphase
} // namespace Physics
} // namespace SF
