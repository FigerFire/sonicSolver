#pragma once

/// @file SF_summary.h
/// @brief Eulerian 方程组一次装配/求解的聚合诊断。

#include <cstdint>

namespace SF::EulerianEulerian {

struct StepSummary {
    int outerCorrectors = 0;
    int pressureCorrections = 0;
    int pressureIterations = 0;
    double pressureResidual = 0.0;
    double maxAlphaSumError = 0.0;
    double massImbalance = 0.0;
    double totalPhaseMass = 0.0;
    double totalPhaseEnthalpy = 0.0;
    double totalWallHeat = 0.0;
    std::int64_t pressureStructureRebuilds = 0;
    std::int64_t pressureLinearSolves = 0;
    int turbulenceIterations = 0;
    double turbulenceResidual = 0.0;
    std::int64_t turbulenceStructureRebuilds = 0;
    std::int64_t turbulenceLinearSolves = 0;
};

} // namespace SF::EulerianEulerian
