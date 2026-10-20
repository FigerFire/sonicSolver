#pragma once

/// @file SF_methodSelection.h
/// @brief 从唯一 IBMRuntimeConfig 构造正交方法选择。

#include "SF_ibmConfig.h"
#include "SF_immersedSystem.h"

namespace SF::IBM::Common {

FDM::ImmersedMethodSelection selectMethod(
    const IBMRuntimeConfig& config);

} // namespace SF::IBM::Common
