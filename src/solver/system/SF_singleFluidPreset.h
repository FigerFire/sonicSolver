#pragma once

/// @file SF_singleFluidPreset.h
/// @brief Built-in density single-fluid equation preset.

#include "SF_buildRequest.h"
#include "SF_equationContribution.h"

namespace SF::System::Preset {

/// @brief Install selected physical balances against explicitly selected solution targets.
///
/// WHAT selects balances; STATE supplies writable metadata and default targets.
/// Default HOW occurrences remain separate from the equation AST and its catalog.
void installSingleFluid(
    SystemCompositionBuilder& system,
    const SingleFluidPresetSpec& spec,
    const EquationCompositionConfig& composition);

} // namespace SF::System::Preset
