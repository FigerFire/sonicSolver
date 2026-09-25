#pragma once

/// @file SF_ibmSystemContribution.h
/// @brief IBM-owned unknown, constraint, equation, and policy contribution.

namespace SF::FDM { struct ImmersedAlgorithmDescriptor; }
#include "core/system/SF_systemContribution.h"

namespace SF::IBM::SystemContribution {

void contribute(
    System::SystemContribution& system,
    const FDM::ImmersedAlgorithmDescriptor& immersed);

} // namespace SF::IBM::SystemContribution
