#pragma once

/// @file SF_frameProvider.h
/// @brief MRF-owned numerical term registrations.

#include "core/system/SF_sourceProvider.h"

namespace SF::Physics::MRF {

std::vector<System::SourceTermProviderDescriptor> termProviders(
    const std::vector<RotatingSetting>& settings);

} // namespace SF::Physics::MRF
