#pragma once

/// @file SF_singleFluidPreset.h
/// @brief Built-in density single-fluid equation preset.

#include "SF_buildRequest.h"
#include "SF_equationContribution.h"

namespace SF::System::Preset {

/// @brief Install rho/rhoU/rhoE and their physical equations.
///
/// This function selects WHAT is solved. It does not select a time recipe,
/// reconstruction, flux, diffusion stencil, or execution order.
void installSingleFluid(
    SystemCompositionBuilder& system,
    const SingleFluidPresetSpec& spec);

} // namespace SF::System::Preset
