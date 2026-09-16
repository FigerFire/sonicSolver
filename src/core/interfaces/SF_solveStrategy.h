#pragma once

/// @file SF_solveStrategy.h
/// @brief 与方程数学内容正交的类型化求解/耦合策略。

#include <string>

namespace SF::FDM {

/// @brief 一个 solve block 的数值耦合方式。
///
/// 名称仅用于 explain/log；runtime 必须只消费 kind，不能从名称推断执行顺序。
enum class SolveStrategyKind {
    ExplicitTimeIntegration,
    SegregatedPredictor,
    PressureCorrection,
    PressureVelocityCoupling,
    ConstraintSolve,
    MonolithicKKT,
    BoundaryClosure,
    AlgebraicUpdate
};

/// @brief 启动阶段冻结的求解策略描述。
struct SolveStrategyDescriptor {
    SolveStrategyKind kind = SolveStrategyKind::AlgebraicUpdate;
    std::string name;
};

inline const char* toString(SolveStrategyKind kind) {
    switch (kind) {
        case SolveStrategyKind::ExplicitTimeIntegration:
            return "explicitTimeIntegration";
        case SolveStrategyKind::SegregatedPredictor:
            return "segregatedPredictor";
        case SolveStrategyKind::PressureCorrection:
            return "pressureCorrection";
        case SolveStrategyKind::PressureVelocityCoupling:
            return "pressureVelocityCoupling";
        case SolveStrategyKind::ConstraintSolve:
            return "constraintSolve";
        case SolveStrategyKind::MonolithicKKT:
            return "monolithicKKT";
        case SolveStrategyKind::BoundaryClosure:
            return "boundaryClosure";
        case SolveStrategyKind::AlgebraicUpdate:
            return "algebraicUpdate";
    }
    return "unknown";
}

} // namespace SF::FDM
