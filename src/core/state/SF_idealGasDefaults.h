#pragma once

/// @file SF_idealGasDefaults.h
/// @brief 尚未绑定 FluidStateModel 时，legacy 单流体路径使用的显式默认值。

namespace SF {

inline constexpr double DefaultIdealGasGamma = 1.4;
inline constexpr double DefaultIdealGasConstant = 287.05;

} // namespace SF
