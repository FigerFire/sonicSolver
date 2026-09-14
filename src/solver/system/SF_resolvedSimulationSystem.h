#pragma once

/// @file SF_resolvedSimulationSystem.h
/// @brief 启动阶段冻结的求解未知量、方程、约束、联立块和执行要求。

#include "SF_configTypes.h"
#include "SF_immersedSystem.h"

#include <string>
#include <vector>

namespace SF::System {

/// @brief 物理状态族；与压力基/密度基数值算法正交。
enum class PhysicsStateKind {
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

/// @brief 一个启动阶段声明的数学未知量。
struct UnknownDescriptor {
    std::string id;
    std::string name;
    VariableLocation location = VariableLocation::EulerianCell;
    int components = 1;
    OwnershipKind ownership = OwnershipKind::EulerianGlobalDof;
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
    std::string strategy;
    std::vector<std::string> equations;
    std::vector<std::string> constraints;
    std::vector<std::string> unknowns;
};

/// @brief 运行所需 backend 能力及其当前可用性。
struct ExecutionRequirement {
    std::string name;
    bool required = false;
    bool available = false;
    std::string detail;
};

/// @brief 一个 case 最终解析出的完整数学系统。
struct ResolvedSimulationSystem {
    FDM::SolverAlgorithm flow = FDM::SolverAlgorithm::DensityBased;
    PhysicsStateKind physics = PhysicsStateKind::SingleFluid;
    std::string timeIntegrator;
    std::vector<UnknownDescriptor> unknowns;
    std::vector<EquationDescriptor> equations;
    std::vector<ConstraintDescriptor> constraints;
    std::vector<SolveBlock> solveBlocks;
    std::vector<std::string> closures;
    std::vector<std::string> boundaries;
    std::vector<ExecutionRequirement> requirements;
    std::string immersedAlgorithm;
    std::string immersedReference;
    std::string immersedSupport;
    std::string immersedRepresentation;
    std::string immersedEnforcement;
    std::string immersedSolid;
    std::string immersedFunctional;
};

/// @brief Builder 所需、已由 application 明确解析的上下文。
struct BuildRequest {
    PhysicsStateKind physics = PhysicsStateKind::SingleFluid;
    std::vector<std::string> phaseNames;
    bool levelSet = false;
    bool turbulence = false;
    std::string turbulenceModel;
    bool parallel = false;
    bool constraintGlobalDofAvailable = false;
    bool distributedLinearSystemAvailable = false;
    const FDM::ImmersedAlgorithmDescriptor* immersed = nullptr;
};

const char* toString(PhysicsStateKind kind);
const char* toString(VariableLocation location);
const char* toString(OwnershipKind ownership);

} // namespace SF::System
