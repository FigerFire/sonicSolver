#pragma once

/// @file SF_planFragment.h
/// @brief Coupling/model 贡献的执行片段 IR。
///
/// 片段只描述 **顺序与重复**，不描述数学：
///   - Leaf 引用 formulation 声明的 `OperationStage`（executable operation
///     才是"存在哪些 operation"的 authority）；
///   - Leaf 不携带 provider/数学语义；
///   - Planner 只负责合并、排序、重复并解析 stage -> OpId。
///
/// 当 formulation 没有为某个 stage 声明 executable operation 时，片段可以
/// 给出 `missingOperation`：这是一个**具名的缺失标记**（用于 Unsupported
/// 报告），不是被 planner 发明出来的数学操作。

#include "core/system/SF_equationIR.h"
#include "core/system/SF_solveProgram.h"

#include <string>
#include <vector>

namespace SF::System {

struct PlanFragmentNode {
    enum class Kind { Sequence, Loop, Leaf };

    Kind kind = Kind::Sequence;
    std::string id;
    std::string name;
    /// @brief Loop 的重复次数。
    int repetitions = 1;
    /// @brief Sequence/Loop 的 plan node kind（Sequence 或 Loop）。
    PlanNodeKind nodeKind = PlanNodeKind::Sequence;
    /// @brief Leaf 的 plan node kind（Assemble/Solve/Correct/Update/Commit…）。
    PlanNodeKind leafKind = PlanNodeKind::Update;
    /// @brief Leaf 引用的 formulation stage。
    OperationStage stage = OperationStage::Prepare;
    /// @brief 已声明 operation 的直接引用，供一个 stage 有多个相操作时使用。
    OpId operation;
    /// @brief formulation 未声明该 stage 时使用的具名缺失标记。
    std::string missingOperation;
    std::string unsupportedReason;
    /// @brief Leaf 关联的 executable equation（仅用于 explain/诊断）。
    std::string equation;
    std::vector<PlanFragmentNode> children;
};

/// @brief 一个可被 planner 合并/排序的执行片段。
struct PlanFragment {
    std::string id;
    std::string name;
    std::string rootId = "PLAN_ROOT";
    std::vector<std::string> consumedPolicies;
    bool includeTimePolicy = true;
    int priority = 20;
    /// @brief 片段根节点的 children（按顺序追加到 plan root）。
    std::vector<PlanFragmentNode> nodes;
};

} // namespace SF::System
