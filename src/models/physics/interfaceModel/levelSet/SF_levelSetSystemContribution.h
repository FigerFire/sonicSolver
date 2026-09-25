#pragma once

/// @file SF_levelSetSystemContribution.h
/// @brief Level-set-owned mathematical system contribution.

#include "core/system/SF_systemContribution.h"

namespace SF::Physics::InterfaceModels::LevelSetContribution {

struct Spec { bool ghostFluid = false; };

void contribute(System::SystemContribution& system, const Spec& spec);

} // namespace SF::Physics::InterfaceModels::LevelSetContribution
