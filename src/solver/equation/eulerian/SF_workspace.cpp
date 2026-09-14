/// @file SF_workspace.cpp
/// @brief 双欧拉方程临时数组与矩阵结构的可复用工作区。

#include "solver/equation/eulerian/SF_workspace.h"
#include "SF_canonicalFace.h"

#include <stdexcept>

namespace SF::EulerianEulerian {

void PhaseFaceFlux::setupLike(const Field& field) {
    const size_t size = 3U * static_cast<size_t>(field.TotalSize());
    volume.assign(size, 0.0);
    mass.assign(size, 0.0);
}

size_t PhaseFaceFlux::index(const Field& field, int direction,
                            int i, int j, int k) {
    return Numerics::CanonicalFace::index(
        field,direction,i,j,k);
}

void PhaseSolverWorkspace::setupLike(
        const Physics::PhaseSystems::PhaseSystem& system) {
    const Field& field = system.geometry();
    const auto& phases = system.phases();
    predictedVelocity.resize(phases.size());
    previousVelocity.resize(phases.size());
    previousPhaseMass.resize(phases.size());
    previousMomentum.resize(phases.size());
    previousPhaseEnthalpy.resize(phases.size());
    momentumDiagonal.resize(phases.size());
    faceFlux.resize(phases.size());
    for (size_t phase = 0; phase < phases.size(); ++phase) {
        momentumDiagonal[phase].setupLike(
            field, "momentumDiagonal." + phases[phase].name);
        previousPhaseMass[phase].setupLike(
            field, "previousPhaseMass." + phases[phase].name);
        previousPhaseEnthalpy[phase].setupLike(
            field, "previousPhaseEnthalpy." + phases[phase].name);
        faceFlux[phase].setupLike(field);
        for (int component = 0; component < 3; ++component) {
            const std::string suffix(1, "xyz"[component]);
            predictedVelocity[phase][(size_t)component].setupLike(
                field, "predictedU" + suffix + "." + phases[phase].name);
            previousVelocity[phase][(size_t)component].setupLike(
                field, "previousU" + suffix + "." + phases[phase].name);
            previousMomentum[phase][(size_t)component].setupLike(
                field, "previousMomentum" + suffix + "."
                + phases[phase].name);
        }
    }
    pressureCorrection.setupLike(field, "pCorrection", 0.0);
    previousPressure.setupLike(field, "previousPressure", 0.0);
    commitTimeLevel(system);
}

void PhaseSolverWorkspace::commitTimeLevel(
        const Physics::PhaseSystems::PhaseSystem& system) {
    if (previousVelocity.size() != system.phases().size()) {
        throw std::runtime_error(
            "PhaseSolverWorkspace is not initialized for this PhaseSystem.");
    }
    for (size_t phase = 0; phase < system.phases().size(); ++phase) {
        previousPhaseMass[phase].values() =
            system.phases()[phase].primary.phaseMass.values();
        previousPhaseEnthalpy[phase].values() =
            system.phases()[phase].primary.phaseEnthalpy.values();
        for (int component = 0; component < 3; ++component) {
            previousVelocity[phase][(size_t)component].values() =
                system.phases()[phase].primitive.velocity[(size_t)component]
                    .values();
            previousMomentum[phase][(size_t)component].values() =
                system.phases()[phase].primary.momentum[(size_t)component]
                    .values();
        }
    }
    previousPressure.values() = system.sharedPressure().values();
}

} // namespace SF::EulerianEulerian
