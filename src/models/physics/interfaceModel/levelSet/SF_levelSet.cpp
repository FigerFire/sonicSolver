/// @file SF_levelSet.cpp
/// @brief Level Set 输运、重初始化或界面几何模型实现。

#include "SF_levelSet.h"

#include "SF_advect.h"
#include "SF_band.h"
#include "SF_heaviside.h"
#include "methods/numerics/structured/SF_structured.h"
#include "SF_surfaceTension.h"
#include "SF_utility.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace SF::Physics::InterfaceModels {
namespace {

double cellVolume(const Field& field, int i, int j, int k) {
    const double inverseJacobian = field.Jac(i, j, k);
    if (!std::isfinite(inverseJacobian)
        || std::abs(inverseJacobian) <= 1.0e-300) {
        std::ostringstream os;
        os << "LevelSet diagnostics: invalid inverse Jacobian at ("
           << i << "," << j << "," << k << ") value="
           << inverseJacobian;
        throw std::runtime_error(os.str());
    }
    return std::abs(1.0 / inverseJacobian);
}

double relativeError(double value, double reference) {
    return value / std::max(std::abs(reference), 1.0e-300);
}

const Multiphase::PhaseProperties* liquidPhase(
        const Multiphase::MultiPhaseConfig& config) {
    for (const auto& phase : config.phases) {
        if (phase.role == Multiphase::PhaseRole::Liquid) return &phase;
    }
    return nullptr;
}

double liquidDensity(const Multiphase::MultiPhaseConfig& config) {
    const auto* phase = liquidPhase(config);
    if (phase == nullptr || !std::isfinite(phase->density)
        || phase->density <= 0.0) {
        throw std::runtime_error(
            "LevelSet diagnostics require one liquid phase with finite "
            "density > 0.");
    }
    return phase->density;
}

double indicator(double phi, double epsilon) {
    return epsilon > 0.0
        ? Multiphase::Heaviside::regularized(phi, epsilon)
        : Multiphase::Heaviside::sharp(phi);
}

} // namespace

LevelSetModel::LevelSetModel(
        const Multiphase::MultiPhaseConfig& config)
    : config_(config) {
    Multiphase::validateMultiPhaseConfig(config_, "LevelSetModel");
    if (!Multiphase::isLevelSetType(config_.type)) {
        throw std::runtime_error(
            "LevelSetModel requires phaseSystem/type = levelSet.");
    }
}

void LevelSetModel::requireInitialized(const char* operation) const {
    if (!initialized_) {
        throw std::runtime_error(
            std::string("LevelSetModel::") + operation
            + " called before initialization.");
    }
}

void LevelSetModel::initialize(const Field& field) {
    state_.initializeFromSets(field, config_);
    state_.applyBoundaryConditions(field, config_);
    state_.computeGeometry(
        field, config_.levelSet.geometryGradientTolerance);
    state_.updateMaterialProperties(field, config_);
    if (config_.levelSet.narrowBandWidth > 0.0) {
        Multiphase::NarrowBand::update(
            field, state_, config_.levelSet.narrowBandWidth);
    }
    rhs_.assign((size_t)state_.TotalSize(), 0.0);
    initialized_ = true;
    updateDiagnostics(field, true);
}

void LevelSetModel::beginTimeStep(double dt) {
    requireInitialized("beginTimeStep");
    if (!std::isfinite(dt) || dt <= 0.0) {
        throw std::runtime_error(
            "LevelSetModel::beginTimeStep requires finite dt > 0.");
    }
    if (config_.levelSet.massCorrection) {
        throw std::runtime_error(
            "LevelSetModel::beginTimeStep: massCorrection is requested, but "
            "global Level Set mass correction is not implemented.");
    }
}

void LevelSetModel::assembleTransportRHS(const Field& field) {
    requireInitialized("assembleTransportRHS");
    Multiphase::AdvectOptions options;
    options.order = config_.levelSet.advectionOrder;
    options.wenoEpsilon = config_.levelSet.wenoEpsilon;
    options.wenoPower = config_.levelSet.wenoPower;
    options.skipSolidCells = config_.levelSet.advectFluidCellsOnly;
    options.skipSolidCellsDeclared = true;
    Multiphase::Advect::assembleRHS(field, state_, options, rhs_);
}

void LevelSetModel::completeTimeStep(
        Field& field, double dt) {
    requireInitialized("completeTimeStep");
    if (!std::isfinite(dt) || dt <= 0.0) {
        throw std::runtime_error(
            "LevelSetModel::completeTimeStep requires finite dt > 0.");
    }
    refresh(field);
}

void LevelSetModel::applyBoundary(const Field& field) {
    requireInitialized("applyBoundary");
    state_.applyBoundaryConditions(field, config_);
}

void LevelSetModel::refresh(const Field& field) {
    requireInitialized("refresh");
    state_.applyBoundaryConditions(field, config_);
    state_.computeGeometry(
        field, config_.levelSet.geometryGradientTolerance);
    state_.updateMaterialProperties(field, config_);
    if (config_.levelSet.narrowBandWidth > 0.0) {
        Multiphase::NarrowBand::update(
            field, state_, config_.levelSet.narrowBandWidth);
    }
    updateDiagnostics(field, false);
}

void LevelSetModel::addSourceTerms(Field& field, Residual& residual) const {
    requireInitialized("addSourceTerms");
    const std::string model = Multiphase::normalizeModelType(
        config_.levelSet.surfaceTensionModel);
    if (model == "none" || model == "ghostfluid") return;
    if (model != "csf") {
        throw std::runtime_error(
            "LevelSetModel::addSourceTerms received an invalid "
            "surfaceTensionModel.");
    }
    if (!state_.hasMaterialProperties()) {
        throw std::runtime_error(
            "LevelSetModel::addSourceTerms requires current material "
            "properties before CSF assembly.");
    }

    const double sigma = config_.levelSet.surfaceTension;
    const double epsilon = config_.levelSet.interfaceThickness;
    Math::forInterior(field, [&](int i, int j, int k) {
        if (field.CellFlag(i, j, k) != FLUID_CELL) return;
        (void)Numerics::requirePhysicalState(
            "LevelSet surface-tension source", field, i, j, k);
        const Vector3 force = Multiphase::SurfaceTension::csfForce(
            state_, i, j, k, sigma, epsilon);
        if (!std::isfinite(force.x) || !std::isfinite(force.y)
            || !std::isfinite(force.z)) {
            std::ostringstream os;
            os << "LevelSet surface tension: non-finite force at ("
               << i << "," << j << "," << k << ").";
            throw std::runtime_error(os.str());
        }
        const double rho = field(i, j, k, RHO);
        const double u = field(i, j, k, RU) / rho;
        const double v = field(i, j, k, RV) / rho;
        const double w = field(i, j, k, RW) / rho;
        residual.source(i, j, k, RU) += force.x;
        residual.source(i, j, k, RV) += force.y;
        residual.source(i, j, k, RW) += force.z;
        residual.source(i, j, k, E) +=
            u * force.x + v * force.y + w * force.z;
    });
}

FDM::InterfacePressureJump LevelSetModel::pressureJump(
        const Field& field, int i, int j, int k, int axis) const {
    requireInitialized("pressureJump");
    FDM::InterfacePressureJump result;
    if (Multiphase::normalizeModelType(
            config_.levelSet.surfaceTensionModel) != "ghostfluid") {
        return result;
    }
    if (!state_.isCompatibleWith(field)) {
        throw std::runtime_error(
            "LevelSetModel::pressureJump requires a compatible Field.");
    }
    if (axis < 0 || axis > 2) {
        throw std::runtime_error(
            "LevelSetModel::pressureJump axis must be 0, 1, or 2.");
    }
    int ri = i, rj = j, rk = k;
    if (axis == 0) ++ri;
    else if (axis == 1) ++rj;
    else ++rk;
    if (i < 0 || i >= field.MX() || j < 0 || j >= field.MY()
        || k < 0 || k >= field.MZ() || ri < 0 || ri >= field.MX()
        || rj < 0 || rj >= field.MY() || rk < 0 || rk >= field.MZ()) {
        throw std::runtime_error(
            "LevelSetModel::pressureJump face index is outside the Field.");
    }

    const double phiLeft = state_.phi(i,j,k);
    const double phiRight = state_.phi(ri,rj,rk);
    if (!std::isfinite(phiLeft) || !std::isfinite(phiRight)) {
        throw std::runtime_error(
            "LevelSetModel::pressureJump found non-finite phi.");
    }
    const bool liquidLeft = phiLeft >= 0.0;
    const bool liquidRight = phiRight >= 0.0;
    if (liquidLeft == liquidRight) return result;

    const double denominator = std::abs(phiLeft) + std::abs(phiRight);
    if (!std::isfinite(denominator) || denominator <= 0.0) {
        throw std::runtime_error(
            "LevelSetModel::pressureJump found an ambiguous zero/zero face.");
    }
    const double fraction = std::abs(phiLeft) / denominator;
    const double curvature =
        (1.0 - fraction) * state_.curvature(i,j,k)
        + fraction * state_.curvature(ri,rj,rk);
    if (!std::isfinite(curvature)) {
        throw std::runtime_error(
            "LevelSetModel::pressureJump found non-finite curvature.");
    }

    // phi>0 为液相，n=grad(phi)/|grad(phi)| 从气相指向液相，
    // 因而 Young-Laplace 约定为 p_liquid-p_gas=-sigma*kappa。
    const double liquidMinusGas =
        -config_.levelSet.surfaceTension * curvature;
    const double phaseOrientation =
        (liquidRight ? 1.0 : 0.0) - (liquidLeft ? 1.0 : 0.0);
    result.crossesInterface = true;
    result.fractionFromLeft = fraction;
    result.targetRightMinusLeft = phaseOrientation * liquidMinusGas;
    if (!std::isfinite(result.targetRightMinusLeft)) {
        throw std::runtime_error(
            "LevelSetModel::pressureJump produced a non-finite jump.");
    }
    return result;
}

void LevelSetModel::correct(const Field& field, double dt) {
    if (!std::isfinite(dt) || dt <= 0.0) {
        throw std::runtime_error(
            "LevelSetModel::correct requires finite dt > 0.");
    }
    refresh(field);
}

double LevelSetModel::dynamicViscosity(
        const Field& field, int i, int j, int k,
        double laminarMu) const {
    (void)laminarMu;
    requireInitialized("dynamicViscosity");
    if (!state_.isCompatibleWith(field)
        || !state_.hasMaterialProperties()) {
        throw std::runtime_error(
            "LevelSetModel::dynamicViscosity requires a current compatible "
            "material-property cache.");
    }
    const double value = state_.viscosity(i, j, k);
    if (!std::isfinite(value) || value < 0.0) {
        std::ostringstream os;
        os << "LevelSetModel::dynamicViscosity: invalid viscosity at ("
           << i << "," << j << "," << k << ") value=" << value;
        throw std::runtime_error(os.str());
    }
    return value;
}

void LevelSetModel::updateDiagnostics(
        const Field& field, bool resetReference) {
    const double rhoLiquid = liquidDensity(config_);
    const double epsilon = config_.levelSet.interfaceThickness;
    double phaseIndicator = 0.0;
    double flowMass = 0.0;
    Math::forInterior(field, [&](int i, int j, int k) {
        if (field.CellFlag(i, j, k) != FLUID_CELL) return;
        const double h = indicator(state_.phi(i, j, k), epsilon);
        const double rho = field(i, j, k, RHO);
        if (!std::isfinite(h) || !std::isfinite(rho) || rho <= 0.0) {
            std::ostringstream os;
            os << "LevelSet diagnostics: invalid state at ("
               << i << "," << j << "," << k << "), H=" << h
               << ", rho=" << rho << ".";
            throw std::runtime_error(os.str());
        }
        const double volume = cellVolume(field, i, j, k);
        phaseIndicator += h * volume;
        flowMass += rho * volume;
    });

    const double phaseMass = phaseIndicator * rhoLiquid;
    if (resetReference || !diagnostics_.initialized) {
        diagnostics_.initialized = true;
        diagnostics_.referencePhaseIndicator = phaseIndicator;
        diagnostics_.referencePhaseMass = phaseMass;
        diagnostics_.referenceFlowMass = flowMass;
    }
    diagnostics_.currentPhaseIndicator = phaseIndicator;
    diagnostics_.phaseIndicatorError =
        phaseIndicator - diagnostics_.referencePhaseIndicator;
    diagnostics_.phaseIndicatorRelativeError = relativeError(
        diagnostics_.phaseIndicatorError,
        diagnostics_.referencePhaseIndicator);
    diagnostics_.currentPhaseMass = phaseMass;
    diagnostics_.phaseMassError = phaseMass - diagnostics_.referencePhaseMass;
    diagnostics_.phaseMassRelativeError = relativeError(
        diagnostics_.phaseMassError, diagnostics_.referencePhaseMass);
    diagnostics_.currentFlowMass = flowMass;
    diagnostics_.flowMassError = flowMass - diagnostics_.referenceFlowMass;
    diagnostics_.flowMassRelativeError = relativeError(
        diagnostics_.flowMassError, diagnostics_.referenceFlowMass);
}

} // namespace SF::Physics::InterfaceModels
