#pragma once

/// @file SF_pressureProviderBinding.h
/// @brief Bind the resolved pressure operation provider once per stepper.

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
    const System::ResolvedOperationBinding& binding,
    const System::CompiledNumericalSystem& numerics,
    const LegacyPressureInputs& legacy);

} // namespace SF::SolverAlgorithm
