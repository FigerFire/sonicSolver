#pragma once

/// @file SF_solveProgram.h
/// @brief ORDER — execution policies、solve blocks 与结构化 control-flow IR。
///        执行顺序只来自 CompiledSolvePlan 的 root。

#include "core/system/SF_equationIR.h"
#include "core/interfaces/SF_solveStrategy.h"

#include <string>
#include <vector>

namespace SF::System {

/// @brief Explain/capability 使用的 typed solve-block view；执行顺序只来自 Plan root。
struct SolveBlock {
    std::string id;
    std::string name;
    std::string strategy;
    std::vector<std::string> equations;
    std::vector<std::string> constraints;
    std::vector<std::string> unknowns;
    FDM::SolveStrategyKind strategyKind = FDM::SolveStrategyKind::AlgebraicUpdate;
};

enum class PlanNodeKind {
    Sequence,
    Loop,
    StageLoop,
    Subcycle,
    BlockSolve,
    Assemble,
    Solve,
    Correct,
    Update,
    Synchronize,
    Reduction,
    Commit,
    ConvergenceCheck
};

/// @brief Open operation identifier resolved by Run::OpRegistry at runtime.
using OpId = std::string;

enum class ExecutionPolicyKind {
    SegregatedPressureCorrection,
    PressureVelocityFixedPoint,
    ConstraintProjection,
    MonolithicKKT,
    BoundaryClosure
};

/// @brief Transformation/preset 对 solve planning 的输入，不含 runtime objects。
struct ExecutionPolicy {
    std::string id;
    std::string name;
    ExecutionPolicyKind kind = ExecutionPolicyKind::BoundaryClosure;
    std::string strategyName;
    FDM::SolveStrategyKind strategyKind =
        FDM::SolveStrategyKind::AlgebraicUpdate;
    std::vector<std::string> equations;
    std::vector<std::string> constraints;
    std::vector<std::string> unknowns;
    int priority = 0;
    int repeatCount = 1;
    int nestedRepeatCount = 1;
    int innerRepeatCount = 1;
    Provenance origin;
    PlanNodeKind leafKind = PlanNodeKind::Update;
    OpId leafOperation;
};

/// @brief Structured control-flow IR；children 的顺序就是执行顺序。
struct SolvePlanNode {
    PlanNodeKind kind = PlanNodeKind::Sequence;
    std::string id;
    std::string name;
    std::vector<std::string> equations;
    std::vector<std::string> constraints;
    OpId operation;
    int repetitions = 1;
    std::vector<SolvePlanNode> children;
    std::string unsupportedReason;
};

/// @brief SolvePlanner 的冻结输出。它只描述执行控制流，不选择 runtime backend。
struct CompiledSolvePlan {
    SolvePlanNode root;
    std::vector<SolveBlock> blocks;
};

/// @brief ORDER/IR 枚举的稳定字符串；header-only，使 run/plan-IR consumers
///        （例如 PlanExecutor）不依赖 system 编译单元。
inline const char* toString(ExecutionPolicyKind value) {
    switch (value) {
        case ExecutionPolicyKind::SegregatedPressureCorrection:
            return "segregated-pressure-correction";
        case ExecutionPolicyKind::PressureVelocityFixedPoint:
            return "pressure-velocity-fixed-point";
        case ExecutionPolicyKind::ConstraintProjection:
            return "constraint-projection";
        case ExecutionPolicyKind::MonolithicKKT: return "monolithic-kkt";
        case ExecutionPolicyKind::BoundaryClosure: return "boundary-closure";
    }
    return "unknown-execution-policy";
}

inline const char* toString(PlanNodeKind value) {
    switch (value) {
        case PlanNodeKind::Sequence: return "Sequence";
        case PlanNodeKind::Loop: return "Loop";
        case PlanNodeKind::StageLoop: return "StageLoop";
        case PlanNodeKind::Subcycle: return "Subcycle";
        case PlanNodeKind::BlockSolve: return "BlockSolve";
        case PlanNodeKind::Assemble: return "Assemble";
        case PlanNodeKind::Solve: return "Solve";
        case PlanNodeKind::Correct: return "Correct";
        case PlanNodeKind::Update: return "Update";
        case PlanNodeKind::Synchronize: return "Synchronize";
        case PlanNodeKind::Reduction: return "Reduction";
        case PlanNodeKind::Commit: return "Commit";
        case PlanNodeKind::ConvergenceCheck: return "ConvergenceCheck";
    }
    return "UnknownPlanNode";
}

} // namespace SF::System
