#pragma once

/// @file SF_capabilityConfig.h
/// @brief 由 IBM 方法选择生成 workflow capability 快照。

#include "SF_immersedSystem.h"

namespace SF::IBM::Common {

FDM::ImmersedMethodCapabilities configureCapabilities(
    const FDM::ImmersedMethodSelection& selection);

} // namespace SF::IBM::Common
