#pragma once

/// @file SF_systemPrinter.h
/// @brief ResolvedSimulationSystem 的稳定启动日志格式化。

#include "SF_resolvedSimulationSystem.h"

#include <string>

namespace SF::System {

std::string describe(const ResolvedSimulationSystem& system);

} // namespace SF::System
