#pragma once

/// @file SF_law.h
/// @brief 与重构格式无关的物理边界定律标识。

namespace SF::Boundary::Law {

/// @brief 物理边界给出的约束；不携带 Linear、MLS 或 ILW 等格式信息。
enum class Kind {
    FixedValue,
    ZeroGradient,
    Symmetry,
    Empty,
    /// @brief 静止无滑移速度壁面 U=0；不决定密度、压力或热边界。
    NoSlip,
    /// @brief 静止平面滑移速度壁面 U_n=0、dU_t/dn=0。
    Slip
};

} // namespace SF::Boundary::Law
