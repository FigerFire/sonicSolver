#pragma once
/// @file SF_pressureConfigTypes.h
/// @brief 压力基耦合算法配置值对象。

#include "core/config/types/SF_linearSolverConfigTypes.h"

namespace SF::FDM {
/// @brief 压力–速度耦合 preset（SIMPLE/PISO/PIMPLE）。
///
/// 它是 coupling preset，不是求解器家族：它不选择物理方程包，只声明
/// pressure-constraint 耦合流程。是否生效由 resolved equation/constraint
/// structure 决定（见 SF_pressureCoupling.h）。
enum class PressureCouplingPreset { SIMPLE, PISO, PIMPLE };

inline const char* toString(PressureCouplingPreset preset) {
    switch (preset) {
        case PressureCouplingPreset::SIMPLE: return "SIMPLE";
        case PressureCouplingPreset::PISO: return "PISO";
        case PressureCouplingPreset::PIMPLE: return "PIMPLE";
    }
    throw std::runtime_error("Unknown pressure coupling algorithm.");
}

/// @brief Eulerian 相级守恒输运的面值格式。
enum class PhaseConvectionScheme { Upwind };

inline const char* toString(PhaseConvectionScheme scheme) {
    if (scheme == PhaseConvectionScheme::Upwind) {
        return "conservativeUpwind";
    }
    throw std::runtime_error("Unknown Eulerian phase convection scheme.");
}

/// @brief 压力基工作流和各方程线性求解配置。
struct PressureCouplingConfig {
    PressureCouplingPreset preset = PressureCouplingPreset::PIMPLE;
    int outerCorrectors = 1;
    int pressureCorrectors = 1;
    int nonOrthogonalCorrectors = 0;
    double momentumRelaxation = 1.0;
    double pressureRelaxation = 1.0;
};

struct PhaseTransportConfig {
    PhaseConvectionScheme convection = PhaseConvectionScheme::Upwind;
    double sourceCfl = 0.5;
};

struct PressureReference {
    int referenceCell = -1;
    double referencePressure = 0.0;
};

struct LinearSolverSet {
    LinearSolverConfig pressure;
    LinearSolverConfig momentum;
    LinearSolverConfig energy;
    LinearSolverConfig turbulence;
};

/// @brief 压力修正控制参数。
struct PressureCorrectionConfig {
    PressureCouplingConfig coupling;
    PhaseTransportConfig phaseTransport;
    PressureReference reference;
    LinearSolverSet linear;
    int maxIterations = 3000;
    double relativeTolerance = 1.0e-2;
    double absoluteTolerance = 1.0e-10;
    double relaxation = 0.001;
    double velocityRelaxation = 0.1;
};

/// @brief 校验 SIMPLE/PISO/PIMPLE 与压力参考配置。
///
/// typed 值对象的自校验；解析层只负责把输入文本变成这些值。
/// @param config 待校验的压力基工作流配置。
inline void validatePressureCorrectionConfig(const PressureCorrectionConfig& config) {
    const auto& coupling = config.coupling;
    if (coupling.outerCorrectors <= 0 || coupling.pressureCorrectors <= 0
        || coupling.nonOrthogonalCorrectors < 0
        || !std::isfinite(coupling.momentumRelaxation)
        || coupling.momentumRelaxation <= 0.0
        || coupling.momentumRelaxation > 1.0
        || !std::isfinite(coupling.pressureRelaxation)
        || coupling.pressureRelaxation <= 0.0
        || coupling.pressureRelaxation > 1.0
        || !std::isfinite(config.phaseTransport.sourceCfl)
        || config.phaseTransport.sourceCfl <= 0.0
        || config.phaseTransport.sourceCfl > 1.0) {
        throw std::invalid_argument(
            "solverProperties has invalid corrector or relaxation controls.");
    }
    if (config.reference.referenceCell < 0
        || !std::isfinite(config.reference.referencePressure)
        || config.reference.referencePressure <= 0.0) {
        throw std::invalid_argument(
            "pressure workflow requires explicit non-negative referenceCell "
            "and finite referencePressure > 0.");
    }
    validateLinearSolverConfig(config.linear.pressure, "pressure linear solver");
    validateLinearSolverConfig(config.linear.momentum, "momentum linear solver");
    validateLinearSolverConfig(config.linear.energy, "energy linear solver");
    validateLinearSolverConfig(config.linear.turbulence,
                               "turbulence linear solver");
}
} // namespace SF::FDM
