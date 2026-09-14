#pragma once

/// @file SF_kkt.h
/// @brief 表面 IBM 的压力—速度—乘子—刚体单体 KKT 校正器。

#include "SF_config.h"
#include "SF_immersedConstraint.h"
#include "solver/linearAlgebra/SF_globalDofSystem.h"

#include <memory>

namespace SF::PressureBased {

/// @brief 一次单体 KKT 求解的收敛与约束诊断。
struct KKTCorrectionSummary {
    int iterations = 0;
    double relativeResidual = 0.0;
    double maxDivergenceBefore = 0.0;
    double maxDivergenceAfter = 0.0;
    std::int64_t structureRebuilds = 0;
    std::int64_t linearSolves = 0;
    FDM::ImmersedConstraintResult immersed;
};

/// @brief 直接装配并求解 `(p',deltaU,Lambda_s,q)` 鞍点系统。
class MonolithicKKT {
public:
    MonolithicKKT(FDM::PressureCorrectionConfig config,
                  double idealGasGamma);

    KKTCorrectionSummary correct(
        Field& field,
        double targetTime,
        double dt,
        FDM::IImmersedConstraint& constraint);

private:
    FDM::PressureCorrectionConfig config_;
    double idealGasGamma_ = 1.4;
    std::vector<LinearAlgebra::GlobalDofId> numberedDofs_;
    std::unique_ptr<LinearAlgebra::StaticDistributedNumbering> numbering_;
    std::unique_ptr<LinearAlgebra::DistributedLinearSystem> linearSystem_;

    /// @brief 为当前数学自由度建立与后端行号隔离的串行编号层。
    void prepareLinearSystem(
        const std::vector<LinearAlgebra::GlobalDofId>& dofs);
};

} // namespace SF::PressureBased
