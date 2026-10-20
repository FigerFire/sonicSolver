#pragma once

/// @file SF_wallFluxProvider.h
/// @brief Wall-heat-owned energy term registration.

#include "core/system/SF_sourceProvider.h"

namespace SF::Physics::WallHeat {

System::SourceTermProviderDescriptor termProvider(
    const std::vector<WallHeatSetting>& settings);

} // namespace SF::Physics::WallHeat
