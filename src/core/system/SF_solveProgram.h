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

using EquationRef = std::string;
using Order = int;

/// Source HOW contains semantics only; providers realize storage downstream.
struct Target {
    std::string symbol;
    TargetKind kind = TargetKind::Physical;
};

inline std::string targetText(const Target& target) {
    return target.symbol + (target.kind==TargetKind::Working ? "*" : "");
}

/// Parse source notation before compilation; kernels receive only typed targets.
inline Target targetFromSyntax(std::string symbol) {
    Target target{std::move(symbol)};
    if(!target.symbol.empty() && target.symbol.back()=='\'')
        throw std::runtime_error("OldTime prime is read-only; correction requires an independently registered delta symbol, not a prime target.");
    if (!target.symbol.empty() && target.symbol.back()=='*') {
        target.kind=TargetKind::Working;
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
    /// @brief HOW compiler根据所属StageLoop冻结；不是第二个时间配置。
    bool stageScoped=false;
    /// @brief 由 STATE 调用绑定冻结，WHICH 不再选择读取版本。
    std::vector<StateUse> inputs;
};

/// @brief 收敛量的数学含义由所选实现提供，控制层只规定范数与判据。
enum class ResidualKind { StateChange, DiscreteEquation, LinearSolver };
enum class ResidualNorm { Maximum, L2 };
enum class ResidualNormalization { None, Initial, Scale };
inline const char* toString(ResidualKind value) {return value==ResidualKind::StateChange?"StateChange":value==ResidualKind::DiscreteEquation?"DiscreteEquation":"LinearSolver";}
inline const char* toString(ResidualNorm value) {return value==ResidualNorm::Maximum?"Maximum":"L2";}
inline const char* toString(ResidualNormalization value) {return value==ResidualNormalization::None?"None":value==ResidualNormalization::Initial?"Initial":"Scale";}
struct ResidualDefinition {
    std::string id,occurrence;
    ResidualKind kind=ResidualKind::StateChange;
    ResidualNorm norm=ResidualNorm::Maximum;
    ResidualNormalization normalization=ResidualNormalization::Initial;
    double absoluteTolerance=1e-10,relativeTolerance=1e-8,scale=1.0;
    /// @brief 编译器冻结采样者，runtime 不按目标名称重新选择实现。
    std::string producerOperation,provider,target;
};

/// @brief HOW 的执行频率；Inherit 继承显式父scope，绝不按provider或相邻节点推断。
enum class ExecutionFrequency { Inherit, Step, Stage };
inline const char* toString(ExecutionFrequency f){return f==ExecutionFrequency::Step?"Step":f==ExecutionFrequency::Stage?"Stage":"Inherit";}

enum class ExecutionKind { Sequence, EquationCall, Loop, StageLoop, Commit, BlockSolve };

struct ExecutionScope {
    ExecutionKind kind = ExecutionKind::Sequence;
    std::string id;
    EquationCall step;
    /// @brief 整体隐式求解的集合；不得逐项zip或展开成Sequence。
    std::vector<std::string> equations,unknowns;
    std::vector<ExecutionScope> children;
    Order order = 0;
    Provenance origin;
    int repetitions = 1;
    int minimumIterations = 1;
    std::string terminationSignal;
    std::vector<std::string> residuals;
    bool residualAll=true,failOnMaximum=true;
    std::vector<Target> commitSources,commitTargets;
    ExecutionFrequency frequency=ExecutionFrequency::Inherit;
};

struct ExecutionProgram {
    /// Root is the single source execution authority.
    ExecutionScope root;
    std::vector<PlacementRequirement> requirements;
    bool explicitOrder=false;
    /// @brief 用户精简 HOW 允许所选 recipe 展开时间 scope；不是 stage 数或时间方法。
    bool recipeScopes=false;
    std::vector<ResidualDefinition> residuals;
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
        case ExecutionKind::BlockSolve: return "BlockSolve";
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
    std::vector<ResidualDefinition> residuals;
    bool residualAll=true,failOnMaximum=true;
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
    /// @brief Ready 后、任何 RHS/代数调用前读取的版本；不自动插入刷新。
    std::vector<StateUse> coefficientInputs;
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

/// @brief 一个已选数值动作的身份；只声明输入准备/结果发布，不携带控制流。
struct BoundOperation { OpId operation;std::string provider; };

struct CompiledEquationCall {
    EquationCall source;
    /// @brief 迁移期间仅用于隔离仍由旧接口编译的调用，不是运行能力或方法选择。
    bool legacyProvider=false;
    CompiledTarget target;
    /// @brief 输入准备按所属物理步/最近迭代/直接非根组解释，不得包含 clock 或树节点。
    std::vector<BoundOperation> stepInputs,iterationInputs,iterationOutputs,groupInputs;
    /// @brief 固定时间 Working predictor 的迭代属于显式信号循环，而非任意内修正循环。
    bool iterationSignalRequired=false;
    /// @brief 本方法的物理结果在 authored Commit 发布；时间推进由过程编译器产生。
    BoundOperation physicalPublication;
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
    /// @brief 显式Stage中的代数求值，不是被时间积分的未知量。
    bool stageLocal=false;
    bool synchronousStages = false;
    std::vector<std::string> stageRequires,stageProvides,coefficientProvides,coefficientStates;
    OpId coefficientPrepare,coefficientPublish;
    StateVersion coefficientVersion=StateVersion::Stage;
    std::vector<StateUse> coefficientInputs;
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
    std::vector<FormulaOperatorBinding> numericalOperators;
    std::vector<Provenance> numericalOrigins;
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
    std::vector<ResidualKind> residualCapabilities;
};

struct CompiledExecutionProgram {
    /// Technical operations compiled from explicit HOW, without another equation.
    std::vector<ExecutableOperation> operations;
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
    std::vector<ResidualDefinition> residuals;
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
