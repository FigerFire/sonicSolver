/// @file SF_prescribed.cpp
/// @brief Bhalla Algorithm 6：给定固体速度的全隐式 DFM/KKT 约束准备。

#include "method/SF_method.h"
#include "SF_config.h"

#include <cmath>
#include <stdexcept>

namespace SF::IBM::Forcing {

const FDM::ImmersedSurfaceSystem&
ImmersedForcingSystem::prepareMonolithicSystem(
        Field& field, double targetTime, double dt) {
    if (!geometry_ || !usesMonolithicKKT()) {
        throw std::runtime_error(
            "ImmersedForcingSystem has no configured fully implicit DLM system.");
    }
    validateImplicitAlgorithmContract();
    if (!std::isfinite(targetTime) || !std::isfinite(dt) || dt <= 0.0) {
        throw std::runtime_error(
            "fullyImplicitDLM requires finite targetTime and positive dt.");
    }
    buildSurfaceSystem(field,targetTime,dt);
    surfaceSystem_.constraintTolerance =
        config_.forcing.constraintTolerance;
    surfaceSystem_.augmentationCoefficient =
        config_.forcing.algorithm
            == FDM::IBMForcingAlgorithm::DFMAugmentedLagrangian
        ? config_.forcing.augmentationCoefficient : 0.0;
    return surfaceSystem_;
}

} // namespace SF::IBM::Forcing
