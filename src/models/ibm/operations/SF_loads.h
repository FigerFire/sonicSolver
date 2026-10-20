#pragma once

/// @file SF_loads.h
/// @brief IBM 乘子到流体/固体载荷及功率的无状态转换。

#include "SF_valueTypes.h"

namespace SF::IBM::Loads {

/// @brief 由流体乘子获得作用在固体上的反力。
Vector3 forceOnBody(const Vector3& fluidMultiplier, double measure);

/// @brief 计算关于质心的固体力矩。
Vector3 torqueOnBody(
    const Vector3& relativePosition,
    const Vector3& bodyForce);

/// @brief 计算乘子对流体所做的机械功率。
double fluidMechanicalPower(
    const Vector3& fluidMultiplier,
    const Vector3& velocity,
    double measure);

} // namespace SF::IBM::Loads
