#pragma once

/// @file SF_resolvedEquationSystem.h
/// @brief WHAT — 主未知量、方程、约束，以及它们到 runtime storage/operator 的
///        编译绑定。Composition → Transformation 的冻结数学 contract。

#include "core/system/SF_expression.h"

#include <string>
#include <string_view>
#include <vector>

namespace SF::System {

enum class VariableLocation {
    EulerianCell,
    EulerianFace,
    BodyConstraint,
    SurfaceConstraint,
    SolidGlobal
};

enum class OwnershipKind {
    EulerianGlobalDof,
    CanonicalFace,
    ConstraintGlobalDof,
    SolidGlobalDof
};

enum class ValueShape { Scalar, Vector, Tensor };
enum class UnknownRole { Primary, Transported, Algebraic, Multiplier, Derived };
enum class StorageBinding {
    PackedDistributed,
    NamedDistributed,
    TransientWorkspace,
    SpecializedExecutor
};

/// @brief Contribution 的来源只服务于 override、validation 与 explain。
enum class OriginKind {
    BuiltinDefault,
    BuiltinPreset,
    Model,
    User,
    Generated
};

/// @brief 用户和 preset 对已有 composition item 的显式修改语义。
enum class ModificationKind { Add, Extend, Replace, Disable };

/// @brief 方程在 raw/executable system 中的数学角色。
enum class EquationCategory {
    PhysicalEquation,
    ConstraintEquation,
    AlgorithmicDerivedEquation,
    AlgebraicRelation
};

struct Provenance {
    OriginKind kind = OriginKind::BuiltinDefault;
    std::string source;
};

struct ContributionRecord {
    std::string id;
    std::string name;
    Provenance origin;
};

struct SystemModification {
    ModificationKind kind = ModificationKind::Add;
    std::string targetKind;
    std::string targetId;
    std::string detail;
    Provenance origin;
};

struct UnknownDescriptor {
    std::string id;
    std::string name;
    VariableLocation location = VariableLocation::EulerianCell;
    int components = 1;
    OwnershipKind ownership = OwnershipKind::EulerianGlobalDof;
    ValueShape shape = ValueShape::Scalar;
    UnknownRole role = UnknownRole::Primary;
    StorageBinding storageBinding = StorageBinding::SpecializedExecutor;
    std::string storageKey;
    int componentOffset = 0;
    bool initializationRequired = true;
    bool boundaryRequired = true;
    bool restartEligible = true;
    bool outputEligible = true;
    bool runtimeStorageRequired = true;
    std::string nameSpace;
    Provenance origin;
};

struct EquationDescriptor {
    std::string id;
    std::string name;
    std::string kind;
    std::vector<std::string> solvedUnknowns;
    EquationCategory category = EquationCategory::PhysicalEquation;
    Provenance origin;
};

struct ConstraintDescriptor {
    std::string id;
    std::string name;
    std::string equation;
    std::string multiplierUnknown;
    Provenance origin;
};

/// @brief Composition 完成、algorithmic transformation 开始前的数学系统。
///
/// 该对象不保存 timestep、RK stage、MPI schedule 或 runner identity。
struct RawEquationSystem {
    std::vector<UnknownDescriptor> unknowns;
    std::vector<EquationDescriptor> equations;
    Equation::System equationDefinitions;
    std::vector<ConstraintDescriptor> constraints;
    std::vector<std::string> closures;
    std::vector<std::string> boundaries;
    std::vector<std::string> dependencies;
    std::vector<ContributionRecord> contributions;
    std::vector<SystemModification> modifications;
};

struct GeneratedOperatorDescriptor {
    std::string id;
    std::string name;
    std::vector<std::string> reads;
    std::vector<std::string> writes;
    Provenance origin;
};

/// @brief Formulation/Transformation 产生的一个 executable operation。
///
/// 它是 "哪些 derived operation 存在" 的唯一 authority：ordering 不在这里，
/// 执行顺序由 solve-plan fragment 决定。runtime 只能执行这里声明过的
/// operation id；provider resolver 负责确认每个 id 都有实现。
enum class OperationStage {
    Prepare,
    MomentumAssemble,
    MomentumSolve,
    PressureBoundaryPrepare,
    PressureAssemble,
    PressureSolve,
    PressureUpdatePrepare,
    VelocityCorrect,
    FluxCorrect,
    CorrectionCommit,
    StepCommit
};

/// @brief Operation 需要的数值能力；不指名具体 implementation provider。
enum class OperationCapability {
    ConservativeExplicit,
    PressureSchedule,
    MomentumPredictor,
    PressureCorrection,
    PressureLinearSolve,
    PressureBoundary,
    VelocityCorrection,
    FluxCorrection,
    EulerianPhaseExecution,
    ImmersedConstraint
};
enum class OperationRecipeRole { Convection, Diffusion };

struct ExecutableOperation {
    /// @brief runtime operation id（provider 必须实现该 id）；值为 `OpId`。
    std::string operation;
    std::string name;
    OperationStage stage = OperationStage::Prepare;
    std::vector<OperationCapability> requirements;
    Provenance origin;
    std::vector<OperationRecipeRole> consumedRecipes;
};

enum class ResourceAccessMode { Read, Write, ReadWrite };
enum class SynchronizationRequirement { None, ReadHalo, WriteOwned };
/// @brief Open execution-provider identifier.  New numerical providers do not
/// require a change to a central enum or to the system compiler.
using OperatorId = std::string;

/// @brief Compiled equation 对 runtime storage/operator 的类型化绑定。
struct CompiledResourceBinding {
    std::string symbol;
    std::string storage;
    int componentOffset = 0;
    int components = 1;
    ResourceAccessMode access = ResourceAccessMode::Read;
    bool boundaryFreshnessRequired = false;
    SynchronizationRequirement synchronization =
        SynchronizationRequirement::None;
};

struct CompiledEquation {
    std::string equationId;
    OperatorId operatorBinding;
    std::vector<CompiledResourceBinding> resources;
    bool assemblesMatrix = false;
    bool assemblesRhs = false;
    Provenance origin;
};

/// @brief Transformation 后、solve planning 前唯一可执行方程 authority。
struct ExecutableEquationSystem {
    std::vector<UnknownDescriptor> unknowns;
    std::vector<EquationDescriptor> equations;
    Equation::System equationDefinitions;
    std::vector<ConstraintDescriptor> constraints;
    /// @brief Formulation 产生的 executable operations（唯一 authority）。
    std::vector<ExecutableOperation> operations;
    std::vector<std::string> closures;
    std::vector<std::string> boundaries;
    std::vector<std::string> dependencies;
    std::vector<GeneratedOperatorDescriptor> correctionOperators;
    std::vector<CompiledEquation> compiledEquations;
};

bool hasUnknown(const RawEquationSystem& system, std::string_view id);
bool hasUnknown(const ExecutableEquationSystem& system, std::string_view id);
bool hasEquation(const RawEquationSystem& system, std::string_view id);
bool hasEquation(const ExecutableEquationSystem& system, std::string_view id);
bool hasConstraint(const RawEquationSystem& system, std::string_view id);
bool hasConstraint(const ExecutableEquationSystem& system, std::string_view id);
bool hasEquationPrefix(const ExecutableEquationSystem& system, std::string_view prefix);

/// @brief 按 stage 查找 executable operation；不存在返回 nullptr。
///
/// Ordering 无关：同一 stage 只允许一个 operation id，plan fragment 只引用
/// stage，具体 id 由 formulation 决定。
const ExecutableOperation* findExecutableOperation(
    const ExecutableEquationSystem& system, OperationStage stage);

const char* toString(OperationStage stage);
const char* toString(OperationCapability capability);

const char* toString(VariableLocation location);
const char* toString(OwnershipKind ownership);
const char* toString(OriginKind kind);
const char* toString(ModificationKind kind);
const char* toString(EquationCategory kind);
const char* toString(ResourceAccessMode kind);
const char* toString(SynchronizationRequirement kind);

} // namespace SF::System
