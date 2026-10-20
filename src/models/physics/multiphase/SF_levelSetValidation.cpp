/// @file SF_levelSetValidation.cpp
/// @brief 多相配置、物性与模型一致性校验实现。

#include "SF_phaseValidation.h"

#include "SF_phaseProperties.h"

#include <cmath>
#include <set>
#include <stdexcept>

namespace SF {
namespace Physics {
namespace Multiphase {

void validateLevelSetConfig(const MultiPhaseConfig& config,
                            const std::string& context) {
    const bool levelSetType = isLevelSetType(config.type);
    if (levelSetType) {
        if (config.defaultPhase.empty()) {
            throw std::runtime_error(
                context + ": [multiPhase].defaultPhase must be explicitly declared.");
        }
        if (!std::isfinite(config.defaultSignedDistance)
            || config.defaultSignedDistance <= 0.0) {
            throw std::runtime_error(
                context + ": [multiPhase].defaultSignedDistance must be finite and > 0.");
        }

        (void)phaseSign(config, config.defaultPhase);

        int liquidCount = 0;
        int gasCount = 0;
        for (const PhaseProperties& phase : config.phases) {
            if (!phase.roleDeclared || phase.role == PhaseRole::Unknown) {
                throw std::runtime_error(
                    context + ": Level Set phase '" + phase.name
                    + "' must explicitly declare role liquid or gas; phase names "
                      "are not used as hidden role defaults.");
            }
            liquidCount += phase.role == PhaseRole::Liquid ? 1 : 0;
            gasCount += phase.role == PhaseRole::Gas ? 1 : 0;
            if (!std::isfinite(phase.density) || phase.density <= 0.0
                || !std::isfinite(phase.viscosity) || phase.viscosity < 0.0) {
                throw std::runtime_error(
                    context + ": Level Set phase '" + phase.name
                    + "' needs explicit SI density > 0 and viscosity >= 0.");
            }
        }
        if (config.phases.size() != 2 || liquidCount != 1 || gasCount != 1) {
            throw std::runtime_error(
                context + ": current one-fluid Level Set requires exactly one "
                  "explicit liquid phase and one explicit gas phase.");
        }

        std::set<std::string> setNames;
        for (const SetPhaseAssignment& assignment : config.setPhases) {
            if (assignment.setName.empty()) {
                throw std::runtime_error(context + ": set name in [sets] cannot be empty.");
            }
            if (assignment.phaseName.empty()) {
                throw std::runtime_error(
                    context + ": set '" + assignment.setName
                    + "' must specify a phase name.");
            }
            if (!setNames.insert(assignment.setName).second) {
                throw std::runtime_error(
                    context + ": duplicate set assignment '" + assignment.setName + "'.");
            }
            (void)phaseSign(config, assignment.phaseName);
        }

        for (const LevelSetPlaneInitializer& plane
             : config.phiPlaneInitializers) {
            if (plane.setName.empty()) {
                throw std::runtime_error(
                    context + ": IC [phi] signed-distance plane set name cannot be empty.");
            }
            const double normalNorm =
                std::sqrt(plane.normal.x * plane.normal.x
                          + plane.normal.y * plane.normal.y
                          + plane.normal.z * plane.normal.z);
            if (!std::isfinite(normalNorm) || normalNorm <= 0.0) {
                throw std::runtime_error(
                    context + ": IC [phi] signed-distance plane set '"
                    + plane.setName + "' needs a finite non-zero normal.");
            }
            if (!std::isfinite(plane.point.x)
                || !std::isfinite(plane.point.y)
                || !std::isfinite(plane.point.z)) {
                throw std::runtime_error(
                    context + ": IC [phi] signed-distance plane set '"
                    + plane.setName + "' point must be finite.");
            }
            if (!std::isfinite(plane.scale) || plane.scale <= 0.0) {
                throw std::runtime_error(
                    context + ": IC [phi] signed-distance plane set '"
                    + plane.setName + "' scale must be finite and > 0.");
            }
        }

        for (const LevelSetSphereInitializer& sphere
             : config.phiSphereInitializers) {
            if (sphere.setName.empty()) {
                throw std::runtime_error(
                    context + ": IC [phi] signed-distance sphere set name cannot be empty.");
            }
            if (!std::isfinite(sphere.center.x)
                || !std::isfinite(sphere.center.y)
                || !std::isfinite(sphere.center.z)
                || !std::isfinite(sphere.radius) || sphere.radius <= 0.0
                || !std::isfinite(sphere.scale) || sphere.scale <= 0.0) {
                throw std::runtime_error(
                    context + ": IC [phi] signed-distance sphere set '"
                    + sphere.setName
                    + "' needs finite center, radius > 0, and scale > 0.");
            }
            if (sphere.insidePhase.empty()) {
                throw std::runtime_error(
                    context + ": IC [phi] signed-distance sphere set '"
                    + sphere.setName + "' must explicitly declare insidePhase.");
            }
            (void)phaseSign(config, sphere.insidePhase);
        }

        for (const LevelSetScalarCondition& condition
             : config.phiInitialConditions) {
            if (condition.setName.empty()) {
                throw std::runtime_error(context + ": IC [phi] set name cannot be empty.");
            }
            if (!condition.typeDeclared || condition.type != FIXED_VALUE) {
                throw std::runtime_error(
                    context + ": IC [phi] set '" + condition.setName
                    + "' must explicitly use FIXED_VALUE.");
            }
            if (!std::isfinite(condition.value)) {
                throw std::runtime_error(
                    context + ": IC [phi] set '" + condition.setName
                    + "' value must be finite.");
            }
        }

        for (const LevelSetScalarCondition& condition
             : config.phiBoundaryConditions) {
            if (condition.setName.empty()) {
                throw std::runtime_error(context + ": BC [phi] set name cannot be empty.");
            }
            if (!condition.typeDeclared
                || (condition.type != FIXED_VALUE
                && condition.type != ZERO_GRADIENT
                && condition.type != EMPTY)) {
                throw std::runtime_error(
                    context + ": BC [phi] set '" + condition.setName
                    + "' supports only FIXED_VALUE, ZERO_GRADIENT, or EMPTY for now.");
            }
            if (condition.type == FIXED_VALUE && !std::isfinite(condition.value)) {
                throw std::runtime_error(
                    context + ": BC [phi] set '" + condition.setName
                    + "' value must be finite.");
            }
        }

        auto validateOrder = [&](int order, const char* name) {
            if (order != 3 && order != 5 && order != 7) {
                throw std::runtime_error(
                    context + ": [levelSet]." + name
                    + " must be 3, 5, or 7.");
            }
        };

        validateOrder(config.levelSet.advectionOrder, "advectionOrder");
        validateOrder(config.levelSet.reinitializationOrder, "reinitializationOrder");

        if (!std::isfinite(config.levelSet.wenoEpsilon)
            || config.levelSet.wenoEpsilon <= 0.0) {
            throw std::runtime_error(
                context + ": [levelSet].wenoEpsilon is required and must be finite > 0.");
        }
        if (!std::isfinite(config.levelSet.wenoPower)
            || config.levelSet.wenoPower <= 0.0) {
            throw std::runtime_error(
                context + ": [levelSet].wenoPower is required and must be finite > 0.");
        }
        if (!config.levelSet.advectFluidCellsOnlyDeclared) {
            throw std::runtime_error(
                context + ": [levelSet].advectFluidCellsOnly must be explicitly set.");
        }

        if (config.levelSet.reinitializationSteps < 0) {
            throw std::runtime_error(
                context + ": [levelSet].reinitializationSteps is required and must be >= 0.");
        }
        if (!std::isfinite(config.levelSet.pseudoTimeStep)
            || config.levelSet.pseudoTimeStep < 0.0) {
            throw std::runtime_error(
                context + ": [levelSet].pseudoTimeStep is required and must be finite >= 0.");
        }
        if (config.levelSet.reinitializationSteps > 0
            && config.levelSet.pseudoTimeStep <= 0.0) {
            throw std::runtime_error(
                context + ": [levelSet].pseudoTimeStep must be finite and > 0 when "
                "reinitializationSteps > 0.");
        }
        if (!std::isfinite(config.levelSet.signSmoothingFactor)
            || config.levelSet.signSmoothingFactor <= 0.0) {
            throw std::runtime_error(
                context + ": [levelSet].signSmoothingFactor is required and must be finite > 0.");
        }
        if (!std::isfinite(config.levelSet.geometryGradientTolerance)
            || config.levelSet.geometryGradientTolerance <= 0.0) {
            throw std::runtime_error(
                context + ": [levelSet].geometryGradientTolerance is required and must be finite > 0.");
        }
        if (config.levelSet.narrowBandWidth < 0.0
            || !std::isfinite(config.levelSet.narrowBandWidth)) {
            throw std::runtime_error(
                context + ": [levelSet].narrowBandWidth must be finite and >= 0.");
        }
        if (config.levelSet.interfaceThickness < 0.0
            || !std::isfinite(config.levelSet.interfaceThickness)) {
            throw std::runtime_error(
                context + ": [levelSet].interfaceThickness must be finite and >= 0.");
        }
        if (config.levelSet.surfaceTension < 0.0
            || !std::isfinite(config.levelSet.surfaceTension)) {
            throw std::runtime_error(
                context + ": [levelSet].surfaceTension must be finite and >= 0.");
        }
        if (!config.levelSet.surfaceTensionModelDeclared) {
            throw std::runtime_error(
                context + ": [levelSet].surfaceTensionModel must be explicitly "
                "set to none, CSF, or ghostFluid.");
        }
        const std::string tensionModel =
            normalizeModelType(config.levelSet.surfaceTensionModel);
        if (tensionModel != "none" && tensionModel != "csf"
            && tensionModel != "ghostfluid") {
            throw std::runtime_error(
                context + ": [levelSet].surfaceTensionModel supports only "
                "none, CSF, or ghostFluid.");
        }
        if (tensionModel == "none" && config.levelSet.surfaceTension != 0.0) {
            throw std::runtime_error(
                context + ": surfaceTensionModel none requires surfaceTension 0.");
        }
        if (tensionModel == "csf"
            && (config.levelSet.surfaceTension <= 0.0
                || config.levelSet.interfaceThickness <= 0.0)) {
            throw std::runtime_error(
                context + ": surfaceTensionModel CSF requires surfaceTension > 0 "
                "and interfaceThickness > 0.");
        }
        if (tensionModel == "ghostfluid"
            && config.levelSet.interfaceThickness != 0.0) {
            throw std::runtime_error(
                context + ": surfaceTensionModel ghostFluid is a sharp jump "
                "method and requires interfaceThickness 0.");
        }
        if (!config.levelSet.massCorrectionDeclared) {
            throw std::runtime_error(
                context + ": [levelSet].massCorrection must be explicitly set.");
        }
        if (config.levelSet.massCorrection) {
            throw std::runtime_error(
                context + ": [levelSet].massCorrection is parsed but not implemented yet.");
        }
    }
}

} // namespace Multiphase
} // namespace Physics
} // namespace SF
