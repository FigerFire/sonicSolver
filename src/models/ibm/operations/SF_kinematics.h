#pragma once

/// @file SF_kinematics.h
/// @brief IBM 固体广义速度与积分点目标速度的无状态运动学操作。

#include "SF_immersedConstraint.h"
#include "SF_ibmConfig.h"

namespace SF::IBM::Forcing {
class IBodyModel;
}

namespace SF::IBM::Kinematics {

/// @brief 当前时间层的刚体平动与转动速度。
struct SolidKinematics {
    Vector3 linearVelocity;
    Vector3 angularVelocity;
};

/// @brief 从固体状态和 KKT 广义速度恢复完整刚体运动学。
SolidKinematics resolve(
    const Forcing::IBodyModel& body,
    const FDM::ImmersedKKTState& state);

/// @brief 计算表面点 `U_s + omega_s x r + u_def`。
Vector3 targetVelocity(
    const SolidKinematics& solid,
    const FDM::ImmersedSurfacePoint& point,
    const Forcing::IBodyModel& body,
    double time);

} // namespace SF::IBM::Kinematics
