#pragma once

/// @file SF_accelerationProvider.h
/// @brief Gravity-owned numerical term registrations.

#include "core/system/SF_sourceProvider.h"

namespace SF::Physics::Gravity {

std::vector<System::SourceTermProviderDescriptor> termProviders(
    const std::vector<ZoneVectorSetting>& settings);

} // namespace SF::Physics::Gravity
