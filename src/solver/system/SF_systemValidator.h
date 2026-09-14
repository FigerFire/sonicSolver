#pragma once

/// @file SF_systemValidator.h
/// @brief ResolvedSimulationSystem 的引用完整性和 backend 能力校验。

#include "SF_resolvedSimulationSystem.h"

namespace SF::System {

void validate(const ResolvedSimulationSystem& system);

} // namespace SF::System
