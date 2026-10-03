#pragma once

/// @file SF_workspace.h
/// @brief Eulerian 求解临时量与每相 canonical face flux。

#include "SF_phaseSystem.h"

#include <cstddef>
#include <vector>

namespace SF::EulerianEulerian {

/// @brief 每相唯一面通量；volume 和 mass 使用相同 canonical face 索引。
struct PhaseFaceFlux {
    std::vector<double> volume;
    std::vector<double> mass;

    void setupLike(const Field& field);
    static size_t index(const Field& field, int direction,
                        int i, int j, int k);
};

/// @brief 求解器临时量；不进入 Physics::PhaseSystem 的物理主状态。
struct PhaseSolverWorkspace {
    std::vector<Physics::PhaseSystems::PhaseVectorField> predictedVelocity;
    std::vector<Physics::PhaseSystems::PhaseVectorField> previousVelocity;
    std::vector<ScalarField> previousPhaseMass;
    std::vector<Physics::PhaseSystems::PhaseVectorField> previousMomentum;
    std::vector<ScalarField> previousPhaseEnthalpy;
    std::vector<ScalarField> momentumDiagonal;
    std::vector<PhaseFaceFlux> faceFlux;
    ScalarField pressureCorrection;
    ScalarField previousPressure;

    void setupLike(const Physics::PhaseSystems::PhaseSystem& system);
    void commitTimeLevel(const Physics::PhaseSystems::PhaseSystem& system);
};

} // namespace SF::EulerianEulerian
