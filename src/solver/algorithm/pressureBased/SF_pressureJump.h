#pragma once

/// @file SF_pressureJump.h
/// @brief 压力校正离散使用的已知锐界面跳跃变换。

#include "core/interfaces/SF_equationCoupling.h"

namespace SF::PressureBased {

/// @brief 一个正方向面上的压力修正跳跃。
struct FaceCorrectionJump {
    bool active = false;
    double fractionFromLeft = 0.0;
    double targetRightMinusLeft = 0.0;
    double correctionRightMinusLeft = 0.0;
};

/// @brief 把目标总压力跳跃换算成压力修正变量的已知跳跃。
///
/// 压力更新为 `p_new=p_old+relaxation*pCorrection`，所以返回值严格满足
/// `pCorrection(right)-pCorrection(left)=(target-current)/relaxation`。
/// 未提供跳跃模型或该面不穿越界面时返回 `active=false`。
FaceCorrectionJump pressureCorrectionJump(
    const FDM::IInterfaceJumpCondition* condition,
    const Field& field,
    int leftI, int leftJ, int leftK,
    int axis,
    double relaxation);

} // namespace SF::PressureBased
