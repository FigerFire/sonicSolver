#pragma once

/// @file SF_solveProgram.h
/// @brief HOW — ordered equation occurrences and generic scope IR; frozen plan lowering.
///        执行顺序只来自 CompiledSolvePlan 的 root。

#include "core/system/SF_equationIR.h"
#include "core/interfaces/SF_solveStrategy.h"
#include "core/config/types/SF_timeRecipe.h"
#include "core/system/SF_numericalBinding.h"
#include "core/system/SF_stateViews.h"

#include <algorithm>
#include <any>
#include <string>
#include <vector>
#include <initializer_list>
#include <optional>

namespace SF::System {

/// @brief Typed execution storage; qualifiers do not prescribe time discretization.
enum class TargetKind { Physical, Working, Correction, Workspace };
using EquationRef = std::string;
using Order = int;

/// Source HOW contains semantics only; providers realize storage downstream.
struct Target {
    std::string symbol;
    TargetKind kind = TargetKind::Physical;
};

inline std::string targetText(const Target& target) {
    return target.symbol + (target.kind==TargetKind::Working ? "*"
        : target.kind==TargetKind::Correction ? "'" : "");
}

/// Parse source notation before compilation; kernels receive only typed targets.
inline Target targetFromSyntax(std::string symbol) {
    Target target{std::move(symbol)};
    if (!target.symbol.empty() && (target.symbol.back()=='*' || target.symbol.back()=='\'')) {
        target.kind=target.symbol.back()=='*' ? TargetKind::Working : TargetKind::Correction;
        target.symbol.pop_back();
    }
    return target;
}

struct CompiledTarget {
    std::string symbol;
    TargetKind kind = TargetKind::Physical;
    std::string workspace;
    VariableLocation workspaceLocation=VariableLocation::EulerianCell;
    OwnershipKind workspaceOwnership=OwnershipKind::EulerianGlobalDof;
    StateViewOwner viewOwner=StateViewOwner::CompilerWorkspace;
    std::vector<CompiledResourceBinding> resources;
    /// @brief 仅显式派生 Update 可写该生产者的输出；普通赋值不能绕过只读契约。
    std::shared_ptr<const StateProducer> producer;
};

/// One occurrence references one mathematical definition and a semantic output.
struct EquationCall {
    EquationRef equation;
    Target target;
    std::string occurrence;
};

enum class ExecutionKind { Sequence, EquationCall, Loop, StageLoop, Commit };

struct ExecutionScope {
    ExecutionKind kind = ExecutionKind::Sequence;
    std::string id;
    EquationCall step;
    std::vector<ExecutionScope> children;
    Order order = 0;
    Provenance origin;
    int repetitions = 1;
    int minimumIterations = 1;
    std::string terminationSignal;
};

struct ExecutionProgram {
    /// Root is the single source execution authority.
    ExecutionScope root;
    std::vector<PlacementRequirement> requirements;
    bool explicitOrder=false;
    ExecutionProgram() = default;
    ExecutionProgram(std::initializer_list<EquationCall> legacySteps)
    {
        for (const auto& call:legacySteps) {
            ExecutionScope node;
            node.kind=ExecutionKind::EquationCall;
            node.step=call;
            root.children.push_back(std::move(node));
        }
    }
};

inline const char* toString(ExecutionKind value) {
    switch (value) {
        case ExecutionKind::Sequence: return "Sequence";
        case ExecutionKind::EquationCall: return "EquationCall";
        case ExecutionKind::Loop: return "Loop";
        case ExecutionKind::StageLoop: return "StageLoop";
        case ExecutionKind::Commit: return "Commit";
    }
    return "UnknownExecutionNode";
}

/// @brief Insertion order breaks ties, independently within every nested scope.
inline void orderExecution(ExecutionScope& scope) {
    std::stable_sort(scope.children.begin(),scope.children.end(),
        [](const ExecutionScope& a,const ExecutionScope& b) { return a.order < b.order; });
    for (auto& child:scope.children) orderExecution(child);
}

/// @brief A compiled reference to WHAT and its explicit HOW output. The
/// selected EquationMethod, never a second mode flag, owns realization.
struct CompiledMathRef {
    EquationRef equation;
    std::string target;
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
using LoopSignalId = std::string;

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
    /// @brief Optional generic signal checked after each complete Loop body.
    LoopSignalId terminationSignal;
    int minimumIterations = 1;
    /// @brief Frozen mathematical calls consumed by this leaf's numerical
    /// provider; an optimized leaf may execute several calls as one kernel.
    std::vector<CompiledMathRef> equationCalls;
    std::optional<CompiledTarget> target;
    std::string occurrence;
    /// Final implementation identity from WHICH; runtime only looks it up.
    std::string provider;
};

/// @brief Provider 声明允许版本；compiler 冻结唯一 coefficient producer/resource。
/// STATE 只定义物理符号；此记录是数值读取契约，不拥有数组或调度。
struct CoefficientUse {
    std::string symbol;
    std::vector<StateVersion> acceptedVersions;
    std::string producer,resource;
    StateVersion version=StateVersion::Current;
};

/// @brief Numerical owner participating in a common recipe, independent of fusion.
struct TemporalParticipant {
    std::string identity, provider;
    std::vector<CompiledMathRef> calls;
    std::vector<std::string> targets;
    std::vector<OpId> preparation;
    OpId stepSize, snapshot, prepareStage, rhs, advance, publish;
    /// @brief 同一 stage RHS 内按 authored HOW 验证的数值资源依赖。
    std::vector<std::string> stageRequires,stageProvides,coefficientProvides,coefficientStates;
    OpId coefficientPrepare,coefficientPublish;
    StateVersion coefficientVersion=StateVersion::Stage;
    std::vector<CoefficientUse> coefficientUses;
};
inline OpId instanceOperation(std::string_view type,std::string_view occurrence) {
    return std::string(type)+"@"+std::to_string(occurrence.size())+":"+std::string(occurrence);
}
namespace TemporalOps {
inline constexpr const char* Provider="execution.temporal";
inline constexpr const char* Dt="time.group.dt";
inline constexpr const char* Open="time.stage.open";
inline constexpr const char* Ready="time.stage.ready";
inline constexpr const char* RhsReady="time.stage.rhsReady";
inline constexpr const char* Close="time.stage.close";
inline constexpr const char* PublishReady="time.group.publishReady";
}

struct CompiledEquationCall {
    EquationCall source;
    CompiledTarget target;
    bool spatialTerms = false;
    bool primitiveSourceRequired = false;
    std::string residualWorkspace;
    int requiredHaloWidth=0;
    /// Numerical storage owned by this provider, described for explain only.
    std::vector<std::string> numericalWorkspaces;
    std::vector<ExecutableOperation> operations;
    std::vector<FDM::TimeRecipeId> temporalCapabilities;
    std::string equationMethod;
    bool temporalResidual = false;
    bool synchronousStages = false;
    std::vector<std::string> stageRequires,stageProvides,coefficientProvides,coefficientStates;
    OpId coefficientPrepare,coefficientPublish;
    StateVersion coefficientVersion=StateVersion::Stage;
    std::vector<CoefficientUse> coefficientUses;
    /// Provider-owned step setup; temporal recipes own stage topology/math.
    std::vector<OpId> temporalPreparation;
    OpId temporalStepSize;
    OpId temporalSnapshot;
    OpId temporalStagePrepare, temporalRhs, temporalAdvance, temporalPublish;
    std::string temporalOwner;
    std::vector<TemporalParticipant> temporalParticipants;
    /// @brief WHICH declares temporal storage; STATE does not select a time backend.
    std::string oldTimeWorkspace;
    std::string stageWorkspace;
    bool publishesStageToPhysicalTarget=false;
    std::vector<CompiledMathRef> calls;
    std::vector<std::string> reads;
    std::vector<std::string> writes;
    std::vector<std::string> sourceMathInputs;
    std::vector<std::string> requirements;
    std::vector<std::string> operatorBindings;
    std::vector<std::string> workspaceRequires;
    std::vector<std::string> workspaceProvides;
    /// @brief 覆写预测器/校正状态使旧资源失效；不插入计算或刷新。
    std::vector<std::string> workspaceInvalidates;
    /// @brief 物理 Commit 前必须已由显式 occurrence 发布的资源。
    std::vector<std::string> commitRequires;
    /// Provider-owned coefficients/correction arrays persist through inner loops of this step.
    std::vector<std::string> stepWorkspaces;
    /// Method-owned numerical micro-topology. A nonempty fragment supersedes
    /// the compatibility backendOperation leaf.
    SolvePlanNode fragment;
    std::string backendOperation;
    std::string backendProvider;
    std::string temporalMethod;
    /// Explicit local fusion contract; a shared backend/provider alone never fuses calls.
    std::string fusionKey;
    std::vector<CompiledMathRef> fusionMembers;
    /// @brief Provider-owned immutable implementation contract; generic HOW lowering never interprets it.
    /// @brief 方法显式读取的限定版本；compiler 冻结到同一 provider-owned view。
    std::vector<CompiledStateView> viewReads;
    std::vector<StateUse> stateUses;
    std::vector<StateEffect> stateEffects;
    std::vector<std::string> capabilities;
    std::vector<CapabilityBinding> boundCapabilities;
    std::vector<CapabilityRequirement> capabilityRequirements;
    /// @brief 所选实现要求的具体运行端口；资源 owner 在 application 绑定，HOW 不分配。
    std::vector<ResourceViewRequirement> physicalResources;
    std::any providerContract;
};

struct CompiledExecutionProgram {
    /// @brief STATE views resolved from HOW targets and WHICH temporal demands.
    std::vector<CompiledStateView> stateViews;
    ExecutionScope root;
    std::vector<CompiledEquationCall> steps;

    /// The selected TemporalMethod, rather than SolvePlanner, owns this
    /// numerical stage topology for a transient equation body.
    SolvePlanNode temporalRoot;
    bool hasTemporalRoot = false;
    std::vector<TemporalParticipant> temporalParticipants;
    /// Provider lifecycle decorates compiled scopes, never source HOW.
    SolvePlanNode loweredRoot;
};

/// @brief SolvePlanner 的冻结输出。它只描述执行控制流，不选择 runtime backend。
struct CompiledSolvePlan {
    SolvePlanNode root;
    ExecutionProgram sourceProgram;
    CompiledExecutionProgram compiledProgram;
};

/// @brief ORDER/IR 枚举的稳定字符串；header-only，使 run/plan-IR consumers
///        （例如 PlanExecutor）不依赖 system 编译单元。
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
