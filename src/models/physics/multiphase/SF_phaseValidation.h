#pragma once

/// @file SF_phaseValidation.h
/// @brief 多相配置校验的内部职责拆分接口。

#include <string>

namespace SF {
namespace Physics {
namespace Multiphase {

struct MultiPhaseConfig;

/// @brief 严格校验 OneFluid Level Set 的显式用户输入。
void validateLevelSetConfig(const MultiPhaseConfig& config,
                            const std::string& context);

} // namespace Multiphase
} // namespace Physics
} // namespace SF
