#pragma once

/// @file SF_immersedSystem.h
/// @brief 与求解算法解耦的统一浸没系统选择接口。
///
/// 该头文件只描述“选择了什么”和“模块能提供什么”，不包含 Field、MPI、
/// HYPRE 或具体离散格式。几何、传递算子、强制策略和固体模型由各自模块
/// 实现；求解器通过小接口组合它们，而不是为每一种组合建立一个新 solver。

#include "core/config/SF_configTypes.h"

#include <string>
#include <vector>

namespace SF::FDM {

/// @brief IBM 未知量在离散系统中的物理位置。
enum class ImmersedVariableLocation {
    EulerianGlobalDof,
    BodyConstraintDof,
    SurfaceConstraintDof,
    SolidGlobalDof
};

/// @brief IBM 未知量的唯一所有权语义。
enum class ImmersedOwnershipKind {
    EulerianOwner,
    ConstraintOwner,
    SolidOwner
};

/// @brief 一个 IBM 算法显式引入的未知量。
struct ImmersedUnknownDescriptor {
    std::string id;
    std::string name;
    ImmersedVariableLocation location =
        ImmersedVariableLocation::EulerianGlobalDof;
    int components = 1;
    ImmersedOwnershipKind ownership =
        ImmersedOwnershipKind::EulerianOwner;
};

/// @brief 一个由乘子施加的 IBM 运动学约束。
struct ImmersedConstraintDescriptor {
    std::string id;
    std::string name;
    std::string equation;
    std::string multiplierUnknown;
};

/// @brief IBM 或固体动力学向统一数学系统贡献的方程。
struct ImmersedEquationDescriptor {
    std::string id;
    std::string name;
    std::string form;
    std::vector<std::string> solvedUnknowns;
};

/// @brief IBM 方程/约束的数值联立方式。
struct ImmersedSolveBlockDescriptor {
    std::string id;
    std::string name;
    std::string strategy;
    std::vector<std::string> equations;
    std::vector<std::string> unknowns;
    std::vector<std::string> constraints;
};

/// @brief IBM 的离散作用量、耗散和约束项声明。
struct ImmersedVariationalDescriptor {
    bool fluidKineticIncrement = false;
    bool viscousDissipation = false;
    bool externalWork = false;
    bool incompressibilityConstraint = false;
    bool immersedNoSlipConstraint = false;
    bool solidKineticIncrement = false;
    bool solidExternalWork = false;
    std::string stationaryFunctional;
};

/// @brief 一个 IBM algorithm 的完整数学身份。
///
/// descriptor 只声明“解什么”，不执行几何、通信或线性求解。所有运行期
/// unknown/constraint 必须先在这里出现，供系统装配和启动日志消费。
struct ImmersedAlgorithmDescriptor {
    std::string id;
    std::string referenceName;
    IBMConstraintSupport support = IBMConstraintSupport::Body;
    IBMRepresentation representation = IBMRepresentation::EulerianMask;
    IBMEnforcement enforcement = IBMEnforcement::GhostCell;
    IBMSolidModel solid = IBMSolidModel::Prescribed;
    bool introducesMultiplier = false;
    bool monolithic = false;
    ImmersedVariationalDescriptor variational;
    std::vector<ImmersedUnknownDescriptor> unknowns;
    std::vector<ImmersedEquationDescriptor> equations;
    std::vector<ImmersedConstraintDescriptor> constraints;
    std::vector<ImmersedSolveBlockDescriptor> solveBlocks;
};

/// @brief 一个被浸没系统约束的流体速度/能量端口。
///
/// `phaseIndex < 0` 表示单流体；非负值表示 Eulerian--Eulerian 的相索引。
/// 端口只保存语义标识，实际 Field 绑定由 equation/algorithm 层完成。
struct ImmersedFluidPort {
    std::string momentumName = "momentum";
    std::string energyName = "energy";
    int phaseIndex = -1;
    bool constrainVelocity = true;
    bool coupleMechanicalWork = true;
};

/// @brief 统一 IBM 的正交选择轴。
struct ImmersedMethodSelection {
    /// @brief IBM family；ghost 与 constraint forcing 不共享同一数学阶段。
    IBMMethod method = IBMMethod::Ghost;
    IBMConstraintSupport support = IBMConstraintSupport::Body;
    IBMRepresentation representation = IBMRepresentation::EulerianMask;
    IBMEnforcement enforcement = IBMEnforcement::GhostCell;
    IBMSolidModel solid = IBMSolidModel::Prescribed;
    std::vector<ImmersedFluidPort> fluidPorts;

    /// @brief 是否包含表面约束。
    bool hasSurfaceConstraint() const {
        return support == IBMConstraintSupport::Surface
            || support == IBMConstraintSupport::SurfaceAndBody;
    }

    /// @brief 是否包含体约束。
    bool hasBodyConstraint() const {
        return support == IBMConstraintSupport::Body
            || support == IBMConstraintSupport::SurfaceAndBody;
    }
};

/// @brief 模块实际提供的能力，不等于用户选择的目标组合。
struct ImmersedMethodCapabilities {
    bool surfaceConstraint = false;
    bool bodyConstraint = false;
    bool diffuseTransfer = false;
    bool eulerianMask = false;
    bool sharpJump = false;
    bool prescribedSolid = false;
    bool coupledRigid = false;
    bool selfPropelledRigid = false;
    bool densityBased = false;
    bool pressureBased = false;
    bool multiPhase = false;
    bool distributed = false;
};

/// @brief 求解器可见的统一浸没系统生命周期。
///
/// 该接口不暴露“何时通信”或“如何组装矩阵”。Algorithm 仍决定数学顺序，
/// ExecutionRuntime 负责 freshness/halo/GlobalDof 契约，具体 enforcement
/// strategy 负责把约束接入相应 EquationSystem。
class IImmersedSystem {
public:
    virtual ~IImmersedSystem() = default;

    /// @brief 返回解析完成的正交选择。
    virtual const ImmersedMethodSelection& methodSelection() const = 0;

    /// @brief 返回实现能力，供 workflow 在分配状态前进行 fail-fast 校验。
    virtual const ImmersedMethodCapabilities& capabilities() const = 0;

    /// @brief 返回在启动阶段冻结的 IBM 数学描述。
    virtual const ImmersedAlgorithmDescriptor& algorithmDescriptor() const = 0;

    /// @brief 返回受约束的流体端口；空列表表示单流体默认端口。
    virtual const std::vector<ImmersedFluidPort>& fluidPorts() const = 0;
};

} // namespace SF::FDM
