#pragma once

/// @file SF_recovery.h
/// @brief 单体 DLM/KKT 解到 Field、乘子场和固体载荷的恢复接口。

#include "SF_ibmConfig.h"
#include "SF_immersedConstraint.h"
#include "operations/SF_kinematics.h"

#include <vector>

namespace SF::IBM::Forcing {
class IBodyModel;
}

namespace SF::IBM::Monolithic {

/// @brief KKT 解恢复产生的完整提交对象。
struct RecoveryResult {
    std::vector<unsigned char> mask;
    std::vector<Vector3> multiplier;
    FDM::ImmersedConstraintResult constraint;
    Kinematics::SolidKinematics solid;
};

/// @brief 恢复表面乘子、载荷、功率和约束残差。
///
/// 该函数不更新刚体时间层；调用者只有在全部校验通过后才能 commit solid。
RecoveryResult recoverSolution(
    Field& field,
    const FDM::ImmersedSurfaceSystem& system,
    const FDM::ImmersedKKTState& state,
    const FDM::IBMForcingConfig& config,
    const Forcing::IBodyModel& body,
    double targetTime);

} // namespace SF::IBM::Monolithic
