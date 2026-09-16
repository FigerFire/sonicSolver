#pragma once

/// @file SF_resolvedSimulationSystem.h
/// @brief 启动阶段冻结的求解未知量、方程、约束、联立块和执行要求。

#include "SF_configTypes.h"
#include "SF_immersedSystem.h"
#include "SF_solveStrategy.h"
#include "solver/equation/SF_expression.h"

#include <string>
#include <string_view>
#include <vector>

namespace SF::System {

/// @brief 物理状态族；与压力基/密度基数值算法正交。
enum class PhysicsTemplateKind {
    SingleFluid,
    HomogeneousMixture,
    OneFluidInterface,
    EulerianEulerian
};

/// @brief 系统未知量的数据位置。
enum class VariableLocation {
    EulerianCell,
    EulerianFace,
    BodyConstraint,
    SurfaceConstraint,
    SolidGlobal
};

/// @brief 系统未知量的 canonical ownership。
enum class OwnershipKind {
    EulerianGlobalDof,
    CanonicalFace,
    ConstraintGlobalDof,
    SolidGlobalDof
};

/// @brief 数学未知量的值形状；独立于底层连续存储布局。
enum class ValueShape { Scalar, Vector, Tensor };

/// @brief 未知量在方程系统中的角色。
enum class UnknownRole { Primary, Transported, Algebraic, Multiplier, Derived };

/// @brief 启动阶段将数学未知量绑定到 runtime storage 的方式。
enum class StorageBinding {
    PackedDistributed,
    NamedDistributed,
    SpecializedExecutor
};

/// @brief 一个启动阶段声明的数学未知量。
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
};

/// @brief 一个独立方程的数学身份。
struct EquationDescriptor {
    std::string id;
    std::string name;
    std::string kind;
    std::vector<std::string> solvedUnknowns;
};

/// @brief 一个乘子约束的数学身份。
struct ConstraintDescriptor {
    std::string id;
    std::string name;
    std::string equation;
    std::string multiplierUnknown;
};

/// @brief 哪些方程、约束和未知量在一个数值阶段共同求解。
struct SolveBlock {
    std::string id;
    std::string name;
    /// @brief 兼容输入/explain 文本；不得被 runtime 用于语义分派。
    std::string strategy;
    std::vector<std::string> equations;
    std::vector<std::string> constraints;
    std::vector<std::string> unknowns;
    /// @brief runtime 只消费的类型化 solve strategy。
    FDM::SolveStrategyKind strategyKind = FDM::SolveStrategyKind::AlgebraicUpdate;
};

/// @brief 运行所需 backend 能力及其当前可用性。
struct ExecutionRequirement {
    std::string name;
    bool required = false;
    bool available = false;
    std::string detail;
};

/// @brief Solver execution 所需、但不属于 physical state 的显式 workspace。
struct WorkspaceRequirement {
    std::string id;
    int components = 1;
    VariableLocation location = VariableLocation::EulerianCell;
    OwnershipKind ownership = OwnershipKind::EulerianGlobalDof;
};

/// @brief 一个 case 最终解析出的完整数学系统。
struct ResolvedSimulationSystem {
    FDM::SolverAlgorithm formulation = FDM::SolverAlgorithm::DensityBased;
    /// @brief 物理状态族仅作为 template provenance / explain 元数据保留。
    ///
    /// Runtime execution 不得以此为 dispatch key 选择 runner；
    /// 执行依据是 unknowns / equations / constraints / solveBlocks。
    PhysicsTemplateKind templateOrigin = PhysicsTemplateKind::SingleFluid;
    std::string timeIntegrator;
    std::vector<UnknownDescriptor> unknowns;
    std::vector<EquationDescriptor> equations;
    /// Runtime assembly 与 check/explain 共享的 executable equation authority。
    Equation::System equationDefinitions;
    std::vector<ConstraintDescriptor> constraints;
    std::vector<SolveBlock> solveBlocks;
    std::vector<std::string> closures;
    std::vector<std::string> boundaries;
    std::vector<ExecutionRequirement> requirements;
    std::vector<WorkspaceRequirement> workspaceRequirements;
    std::string immersedAlgorithm;
    std::string immersedReference;
    std::string immersedSupport;
    std::string immersedRepresentation;
    std::string immersedEnforcement;
    std::string immersedSolid;
    std::string immersedFunctional;
};

/// @brief 只读查询；运行装配应查询数学系统，而不是 template provenance。
bool hasUnknown(const ResolvedSimulationSystem& system, std::string_view id);
bool hasEquation(const ResolvedSimulationSystem& system, std::string_view id);
const Equation::Definition& equationDefinition(
    const ResolvedSimulationSystem& system, std::string_view id);
bool hasEquationPrefix(
    const ResolvedSimulationSystem& system, std::string_view prefix);
bool hasConstraint(const ResolvedSimulationSystem& system, std::string_view id);
bool hasSolveBlock(const ResolvedSimulationSystem& system, std::string_view id);
bool hasSolveStrategy(
    const ResolvedSimulationSystem& system, FDM::SolveStrategyKind kind);
bool hasRequirement(
    const ResolvedSimulationSystem& system, std::string_view name);
bool requiresCapability(
    const ResolvedSimulationSystem& system, std::string_view name);

/// @brief Builder 所需、已由 application 明确解析的上下文。
struct BuildRequest {
    PhysicsTemplateKind templateOrigin = PhysicsTemplateKind::SingleFluid;
    std::vector<std::string> phaseNames;
    bool levelSet = false;
    bool homogeneousThermodynamics = false;
    bool legacyMixture = false;
    bool phaseChange = false;
    bool transportedLegacyAlpha = false;
    bool interfaceGhostFluid = false;
    bool turbulence = false;
    std::string turbulenceModel;
    std::vector<std::string> turbulencePhaseNames;
    bool parallel = false;
    bool constraintGlobalDofAvailable = false;
    bool distributedLinearSystemAvailable = false;
    const FDM::ImmersedAlgorithmDescriptor* immersed = nullptr;
};

const char* toString(PhysicsTemplateKind kind);
const char* toString(VariableLocation location);
const char* toString(OwnershipKind ownership);

} // namespace SF::System
