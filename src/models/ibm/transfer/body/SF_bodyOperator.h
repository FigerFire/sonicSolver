#pragma once

/// @file SF_bodyOperator.h
/// @brief Eulerian body-mask 约束的 canonical 单 patch 拓扑构造器。

#include "SF_field.h"
#include "constraint/SF_constraintOperator.h"
#include "solid/SF_bodyModel.h"

namespace SF::IBM::Forcing {

/// @brief 只负责识别 body constraint DOF 与其离散对偶体积。
///
/// 流体密度、动量、能量、目标速度和乘子均不存入该对象。这样
/// fractional DLM 与 Brinkman enforcement 可以消费完全相同的几何拓扑。
class BodyOperator final : public IBodyConstraintOperator {
public:
    /// @brief 为当前刚体构型构造确定顺序的 Eulerian body constraint set。
    const FDM::ImmersedBodySystem& build(
        const Field& field,
        const GeoProcessing::STLGeometry& geometry,
        const IBodyModel& bodyModel,
        double targetTime) override;

    const FDM::ImmersedBodySystem& system() const { return system_; }

private:
    FDM::ImmersedBodySystem system_;
};

} // namespace SF::IBM::Forcing
