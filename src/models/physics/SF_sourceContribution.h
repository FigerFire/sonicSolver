#pragma once

/// @file SF_sourceContribution.h
/// @brief Physical source models contribute terms to existing equations.

#include "core/config/SF_configTypes.h"

#include "core/system/SF_systemContribution.h"

namespace SF::Physics::SourceContribution {

void contribute(
    System::SystemContribution& system,
    const FDM::SourceConfig& config);

} // namespace SF::Physics::SourceContribution
