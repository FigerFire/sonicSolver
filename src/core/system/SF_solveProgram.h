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
    /// Root is the single source execution authority; legacy entries are explicit migration debt.
    ExecutionScope root;
    std::vector<PlacementRequirement> requirements;
    bool explicitOrder=false;
    std::vector<ExecutionScope> legacyEntries;
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

enum class LegacyExecutionPolicyKind {
    SegregatedPressureCorrection,
    PressureVelocityFixedPoint,
    BoundaryClosure
};

/// @brief Explain/capability 使用的 typed solve-block view；执行顺序只来自 Plan root。
struct SolveBlock {
    std::string id;
    std::string name;
    std::string strategy;
    std::vector<std::string> equations;
    std::vector<std::string> constraints;
    std::vector<std::string> unknowns;
    FDM::SolveStrategyKind strategyKind = FDM::SolveStrategyKind::AlgebraicUpdate;
    LegacyExecutionPolicyKind policyKind = LegacyExecutionPolicyKind::BoundaryClosure;
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

/// @brief Transformation/preset 对 solve planning 的输入，不含 runtime objects。
struct LegacyExecutionPolicy {
    std::string id;
    std::string name;
    LegacyExecutionPolicyKind kind = LegacyExecutionPolicyKind::BoundaryClosure;
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
    /// Only explicitly lowered compatibility leaves may use the legacy adapter.
    bool legacyAdapter = false;
};

/// @brief Numerical owner participating in a common recipe, independent of fusion.
struct TemporalParticipant {
    std::string identity, provider;
    std::vector<CompiledMathRef> calls;
    std::vector<std::string> targets;
    std::vector<OpId> preparation;
    OpId stepSize, snapshot, prepareStage, rhs, advance, publish;
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
    std::vector<StateUse> stateUses;
    std::vector<StateEffect> stateEffects;
    std::vector<std::string> capabilities;
    std::vector<CapabilityBinding> boundCapabilities;
    std::vector<CapabilityRequirement> capabilityRequirements;
    std::any providerContract;
};

struct CompiledExecutionProgram {
    /// @brief STATE views resolved from HOW targets and WHICH temporal demands.
    std::vector<CompiledStateView> stateViews;
    ExecutionScope root;
    std::vector<CompiledEquationCall> steps;
    /// Explicit legacy numerical inputs; never inferred by the generic compiler.
    std::vector<CompiledMathRef> legacySpatialInputs;
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
    std::vector<SolveBlock> blocks;
    ExecutionProgram sourceProgram;
    CompiledExecutionProgram compiledProgram;
};

/// @brief ORDER/IR 枚举的稳定字符串；header-only，使 run/plan-IR consumers
///        （例如 PlanExecutor）不依赖 system 编译单元。
inline const char* toString(LegacyExecutionPolicyKind value) {
    switch (value) {
        case LegacyExecutionPolicyKind::SegregatedPressureCorrection:
            return "segregated-pressure-correction";
        case LegacyExecutionPolicyKind::PressureVelocityFixedPoint:
            return "pressure-velocity-fixed-point";
        case LegacyExecutionPolicyKind::BoundaryClosure: return "boundary-closure";
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
