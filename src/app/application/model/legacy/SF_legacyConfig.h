#pragma once

/// @file SF_legacyConfig.h
/// @brief 将旧 parser 暂存量一次性转换成求解器值对象。

#include "SF_configTypes.h"

namespace SF::Legacy {

/// @brief 复制 parser 全局暂存值，生成不依赖全局状态的求解器配置。
FDM::SolverConfig makeSolverConfig();

} // namespace SF::Legacy
