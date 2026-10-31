#pragma once

/// @file SF_pressureProviderBinding.h
/// @brief 从冻结方法/资源构造压力数值 workspace；conservative corrector 保持独立窄输入。

#include "solver/algorithm/pressure/SF_pressureOperators.h"
#include "solver/algorithm/pressureBased/SF_corrector.h"
#include "solver/system/SF_runtimeRequirements.h"

#include <memory>

namespace SF::SolverAlgorithm {

struct PressureProviderBinding {
    std::unique_ptr<Pressure::PressureOperators> constant;
    std::unique_ptr<PressureBased::Corrector> conservative;
};

/// @brief Transitional inputs used only by the existing conservative corrector.
struct LegacyPressureInputs {
    const FDM::BoundaryConfig& boundaries;
    const FDM::PressureCorrectionConfig& pressure;
    double idealGasGamma;
};

PressureProviderBinding bindPressureProvider(
    const System::Program& program,
    const LegacyPressureInputs& legacy);

} // namespace SF::SolverAlgorithm
