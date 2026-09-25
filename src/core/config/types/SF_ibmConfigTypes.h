#pragma once
/// @file SF_ibmConfigTypes.h
/// @brief IBM 方法、几何分类与 forcing 配置值对象。

#include "core/state/SF_valueTypes.h"

#include <limits>
#include <string>
#include <vector>

namespace SF::FDM {
/// @brief WENO stencil 接触 IBM 时采用的闭合方式。
enum class IBMBoundaryScheme { LowOrder, ILW };
/// @brief IBM 几何分类和 ILW 法向搜索配置。
enum class IBMMethod { Ghost, VariationalForcing };

/// @brief Bhalla 统一约束表述下的可执行 IBM 算法。
///
/// FTS 与 fully implicit 是互斥的时间耦合方式；枚举名称按论文
/// Algorithms 1--7 区分，不把二者合并成含混的“implicit FTS”。
enum class IBMForcingAlgorithm {
    PeskinOriginal,
    DFMExplicitSelfPropelled,
    DFMFractionalStepSelfPropelled,
    DFMFractionalStepPrescribed,
    DFMImplicitPrescribed,
    DFMImplicitSelfPropelled,
    DFMAugmentedLagrangian,
    VelocityForcingFTS,
    VelocityForcingBP
};

/// @brief 变分 IBM 施加约束的几何域。
enum class IBMConstraintDomain { Volume, Surface };

/// @brief 约束支持域。与离散强制方式正交，允许同时声明表面和体约束。
///
/// `IBMConstraintDomain` 保留为旧版输入兼容字段；新接口应优先使用
/// `constraintSupport`，因为 Bhalla 统一表述中的 Lambda_s 与 lambda_b
/// 不是互相替代的两种算法。
enum class IBMConstraintSupport {
    Body,
    Surface,
    SurfaceAndBody
};

/// @brief 浸没力/乘子的离散表示。
enum class IBMRepresentation {
    DiffuseKernel,
    EulerianMask,
    SharpJump
};

/// @brief 浸没约束的时间/代数强制策略。
enum class IBMEnforcement {
    GhostCell,
    ExplicitIBM,
    FractionalDLM,
    VelocityForcing,
    BrinkmanPenalty,
    MonolithicKKT
};

/// @brief 固体未知量的提供方式。
enum class IBMSolidModel {
    Prescribed,
    CoupledRigid,
    SelfPropelledRigid,
    Deformable
};

/// @brief 流体、压力和乘子之间的时间层耦合方式。
enum class IBMConstraintCoupling {
    IncrementalProjection,
    MonolithicKKT
};

/// @brief 固体运动自由度的提供方式。
enum class IBMSolidMotion { PrescribedRigid, CoupledRigid };

/// @brief 刚体方程中启用的运动自由度；motivation 表示平动。
enum class IBMRigidMotionMode { Motivation, Rotate };

/// @brief 约束校正对流体能量的离散处理。
enum class IBMEnergyCoupling { MechanicalWork };

/// @brief 变分/DLM forcing 的显式用户配置。
struct IBMForcingConfig {
    /// @brief 真正选择求解流程的算法键；其余枚举描述空间离散和耦合轴。
    IBMForcingAlgorithm algorithm =
        IBMForcingAlgorithm::DFMFractionalStepPrescribed;
    /// @brief 旧版二选一输入；使用 constraintSupport 后只作兼容映射。
    IBMConstraintDomain constraintDomain = IBMConstraintDomain::Volume;
    /// @brief 新版正交约束域选择。
    IBMConstraintSupport constraintSupport = IBMConstraintSupport::Body;
    /// @brief 约束在 Eulerian 网格上的表示。
    IBMRepresentation representation = IBMRepresentation::EulerianMask;
    /// @brief 乘子/强制的时间耦合策略。
    IBMEnforcement enforcement = IBMEnforcement::FractionalDLM;
    /// @brief 固体自由度模型。
    IBMSolidModel solidModel = IBMSolidModel::Prescribed;
    /// @brief 被约束的流体端口；空值表示单流体 momentum/energy。
    std::vector<std::string> fluidPorts;
    IBMConstraintCoupling coupling =
        IBMConstraintCoupling::IncrementalProjection;
    IBMSolidMotion motion = IBMSolidMotion::PrescribedRigid;
    IBMEnergyCoupling energyCoupling =
        IBMEnergyCoupling::MechanicalWork;
    Vector3 centerOfMass;
    Vector3 linearVelocity;
    Vector3 angularVelocity;
    double constraintTolerance = 1.0e-10;
    /// @brief 串行表面 Schur 投影的最大 CG 迭代数；不提供内置猜测值。
    int constraintSolverMaxIterations = 0;
    /// @brief 串行表面 Schur 投影的相对残差容差。
    double constraintSolverRelativeTolerance =
        std::numeric_limits<double>::quiet_NaN();
    IBMRigidMotionMode rigidMotionMode =
        IBMRigidMotionMode::Motivation;
    std::string surfaceKernel;
    std::string surfaceQuadrature;
    std::string surfaceNormalization;
    std::string surfaceSpreading;
    double surfaceSupportRadius =
        std::numeric_limits<double>::quiet_NaN();
    /// @brief Brinkman penalty 系数；仅 penalty 策略允许读取，单位按方程定义。
    double penaltyCoefficient = std::numeric_limits<double>::quiet_NaN();
    /// @brief Augmented-Lagrangian 的 gamma，单位与离散动量 Hessian 一致。
    double augmentationCoefficient =
        std::numeric_limits<double>::quiet_NaN();
    double rigidMass = std::numeric_limits<double>::quiet_NaN();
    Vector3 principalInertia{
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN()};
    Vector3 externalForce;
    Vector3 externalTorque;
};

/// @brief IBM 方法、几何分类和 ILW/forcing 配置。
struct IBMConfig {
    bool enabled = false;
    IBMMethod method = IBMMethod::Ghost;
    double normalAngleDegrees = 69.51268488527785;
    int normalSearchMinLayers = 1;
    int normalSearchMaxLayers = 0;
    int normalSearchTargetCandidates = 512;
    int normalSearchKeepSamples = 256;
    IBMForcingConfig forcing;
};
inline const char* toString(IBMMethod method) {
    return method == IBMMethod::Ghost ? "ghost" : "variationalForcing";
}
/// @brief 将 Bhalla IBM 算法转换为稳定输入/日志 token。
inline const char* toString(IBMForcingAlgorithm algorithm) {
    switch (algorithm) {
        case IBMForcingAlgorithm::PeskinOriginal:
            return "peskinOriginal";
        case IBMForcingAlgorithm::DFMExplicitSelfPropelled:
            return "dfmExplicitSelfPropelled";
        case IBMForcingAlgorithm::DFMFractionalStepSelfPropelled:
            return "dfmFractionalStepSelfPropelled";
        case IBMForcingAlgorithm::DFMFractionalStepPrescribed:
            return "dfmFractionalStepPrescribed";
        case IBMForcingAlgorithm::DFMImplicitPrescribed:
            return "dfmImplicitPrescribed";
        case IBMForcingAlgorithm::DFMImplicitSelfPropelled:
            return "dfmImplicitSelfPropelled";
        case IBMForcingAlgorithm::DFMAugmentedLagrangian:
            return "dfmAugmentedLagrangian";
        case IBMForcingAlgorithm::VelocityForcingFTS:
            return "velocityForcingFTS";
        case IBMForcingAlgorithm::VelocityForcingBP:
            return "velocityForcingBP";
    }
    return "dfmFractionalStepPrescribed";
}
/// @brief 将浸没约束支持域转换为稳定日志 token。
inline const char* toString(IBMConstraintSupport support) {
    switch (support) {
        case IBMConstraintSupport::Body: return "body";
        case IBMConstraintSupport::Surface: return "surface";
        case IBMConstraintSupport::SurfaceAndBody: return "surfaceAndBody";
    }
    return "body";
}
/// @brief 将浸没力空间表示转换为稳定日志 token。
inline const char* toString(IBMRepresentation representation) {
    switch (representation) {
        case IBMRepresentation::DiffuseKernel: return "diffuseKernel";
        case IBMRepresentation::EulerianMask: return "eulerianMask";
        case IBMRepresentation::SharpJump: return "sharpJump";
    }
    return "diffuseKernel";
}
/// @brief 将浸没强制策略转换为稳定日志 token。
inline const char* toString(IBMEnforcement enforcement) {
    switch (enforcement) {
        case IBMEnforcement::GhostCell: return "ghostCell";
        case IBMEnforcement::ExplicitIBM: return "explicitIBM";
        case IBMEnforcement::FractionalDLM: return "fractionalDLM";
        case IBMEnforcement::VelocityForcing: return "velocityForcing";
        case IBMEnforcement::BrinkmanPenalty: return "brinkmanPenalty";
        case IBMEnforcement::MonolithicKKT: return "monolithicKKT";
    }
    return "fractionalDLM";
}
/// @brief 将固体模型转换为稳定日志 token。
inline const char* toString(IBMSolidModel model) {
    switch (model) {
        case IBMSolidModel::Prescribed: return "prescribed";
        case IBMSolidModel::CoupledRigid: return "coupledRigid";
        case IBMSolidModel::SelfPropelledRigid: return "selfPropelledRigid";
        case IBMSolidModel::Deformable: return "deformable";
    }
    return "prescribed";
}
/// @brief Convert an IBM boundary closure mode to the canonical config/log token.
/// @param scheme Strongly typed IBM boundary scheme.
/// @return Stable string used in logs and dispatch bridges.
inline const char* toString(IBMBoundaryScheme scheme) {
    switch (scheme) {
        case IBMBoundaryScheme::LowOrder: return "Downgrade";
        case IBMBoundaryScheme::ILW: return "ILW";
    }
    return "Downgrade";
}
inline void validateIBMForcingAlgorithm(const IBMForcingConfig& config) {
    const bool body = config.constraintSupport == IBMConstraintSupport::Body;
    const bool surface =
        config.constraintSupport == IBMConstraintSupport::Surface;
    const auto requireBodyFTS = [&]() {
        if (!body || config.representation != IBMRepresentation::EulerianMask
            || config.enforcement != IBMEnforcement::FractionalDLM
            || config.coupling
                != IBMConstraintCoupling::IncrementalProjection) {
            throw std::invalid_argument(
                "DFM fractional-step algorithms require body + "
                "eulerianMask + fractionalDLM + incrementalProjection.");
        }
    };
    const auto requireSurfaceKKT = [&]() {
        if (!surface
            || config.representation != IBMRepresentation::DiffuseKernel
            || config.enforcement != IBMEnforcement::MonolithicKKT
            || config.coupling != IBMConstraintCoupling::MonolithicKKT) {
            throw std::invalid_argument(
                "DFM implicit algorithms require surface + diffuseKernel + "
                "monolithicKKT coupling.");
        }
    };
    switch (config.algorithm) {
        case IBMForcingAlgorithm::PeskinOriginal:
            if (!surface
                || config.representation != IBMRepresentation::DiffuseKernel
                || config.enforcement != IBMEnforcement::ExplicitIBM
                || config.solidModel != IBMSolidModel::Prescribed) {
                throw std::invalid_argument(
                    "peskinOriginal requires surface + diffuseKernel + "
                    "explicitIBM + prescribed solidModel. The elastic-stress "
                    "provider required by a deformable Peskin solid is not "
                    "implemented yet.");
            }
            break;
        case IBMForcingAlgorithm::DFMExplicitSelfPropelled:
            if (!body || config.representation != IBMRepresentation::EulerianMask
                || config.enforcement != IBMEnforcement::ExplicitIBM
                || config.solidModel != IBMSolidModel::SelfPropelledRigid) {
                throw std::invalid_argument(
                    "dfmExplicitSelfPropelled requires body + eulerianMask + "
                    "explicitIBM + selfPropelledRigid.");
            }
            break;
        case IBMForcingAlgorithm::DFMFractionalStepSelfPropelled:
            requireBodyFTS();
            if (config.solidModel != IBMSolidModel::SelfPropelledRigid) {
                throw std::invalid_argument(
                    "dfmFractionalStepSelfPropelled requires "
                    "solidModel selfPropelledRigid.");
            }
            break;
        case IBMForcingAlgorithm::DFMFractionalStepPrescribed:
            requireBodyFTS();
            if (config.solidModel != IBMSolidModel::Prescribed) {
                throw std::invalid_argument(
                    "dfmFractionalStepPrescribed requires solidModel prescribed.");
            }
            break;
        case IBMForcingAlgorithm::DFMImplicitPrescribed:
            requireSurfaceKKT();
            if (config.solidModel != IBMSolidModel::Prescribed
                && config.solidModel != IBMSolidModel::CoupledRigid) {
                throw std::invalid_argument(
                    "dfmImplicitPrescribed requires prescribed or "
                    "coupledRigid solidModel.");
            }
            break;
        case IBMForcingAlgorithm::DFMImplicitSelfPropelled:
            requireSurfaceKKT();
            if (config.solidModel != IBMSolidModel::SelfPropelledRigid) {
                throw std::invalid_argument(
                    "dfmImplicitSelfPropelled requires "
                    "solidModel selfPropelledRigid.");
            }
            break;
        case IBMForcingAlgorithm::DFMAugmentedLagrangian:
            requireSurfaceKKT();
            if (config.solidModel != IBMSolidModel::Prescribed) {
                throw std::invalid_argument(
                    "dfmAugmentedLagrangian currently requires prescribed "
                    "solidModel; coupled-solid augmentation needs the full "
                    "u-q cross block and is not silently approximated.");
            }
            if (!std::isfinite(config.augmentationCoefficient)
                || config.augmentationCoefficient <= 0.0) {
                throw std::invalid_argument(
                    "dfmAugmentedLagrangian requires a finite positive "
                    "augmentationCoefficient.");
            }
            break;
        case IBMForcingAlgorithm::VelocityForcingFTS:
            if (!surface
                || config.representation != IBMRepresentation::DiffuseKernel
                || config.enforcement != IBMEnforcement::VelocityForcing) {
                throw std::invalid_argument(
                    "velocityForcingFTS requires surface + diffuseKernel + "
                    "velocityForcing.");
            }
            if (config.constraintSolverMaxIterations <= 0
                || !std::isfinite(
                    config.constraintSolverRelativeTolerance)
                || config.constraintSolverRelativeTolerance <= 0.0) {
                throw std::invalid_argument(
                    "velocityForcingFTS requires explicit positive "
                    "constraintSolver maxIterations and relativeTolerance.");
            }
            break;
        case IBMForcingAlgorithm::VelocityForcingBP:
            if (!body || config.representation != IBMRepresentation::EulerianMask
                || config.enforcement != IBMEnforcement::BrinkmanPenalty) {
                throw std::invalid_argument(
                    "velocityForcingBP requires body + eulerianMask + "
                    "brinkmanPenalty.");
            }
            break;
    }
}

/// @brief 校验 algorithm 与空间/时间离散轴是否描述同一方法。
// 校验已移入 core/config/types/SF_ibmConfigTypes.h：algorithm/selection 校验是
// IBMForcingConfig 的 typed 值对象自校验，models/solver 不能反向依赖
// application 解析层。
/// @brief 校验 IBM 正交选择的物理与实现边界。
///
/// 该函数只检查选择之间的显式一致性；它不会把一种策略替换成另一种
/// 策略。尚未迁移到求解器的组合在构造 forcing 时继续给出专门诊断。
inline void validateIBMForcingSelection(const IBMForcingConfig& config) {
    if (config.enforcement == IBMEnforcement::GhostCell) {
        throw std::invalid_argument(
            "ghostCell is an independent IBM method and cannot be selected "
            "inside a variational forcing block.");
    }
    const bool surface =
        config.constraintSupport == IBMConstraintSupport::Surface
        || config.constraintSupport == IBMConstraintSupport::SurfaceAndBody;
    const bool body =
        config.constraintSupport == IBMConstraintSupport::Body
        || config.constraintSupport == IBMConstraintSupport::SurfaceAndBody;

    if (!surface && !body) {
        throw std::invalid_argument(
            "IBM forcing must select at least one constraint support domain.");
    }
    if (config.representation == IBMRepresentation::EulerianMask && !body) {
        throw std::invalid_argument(
            "eulerianMask representation requires a body/volume constraint.");
    }
    if (config.representation == IBMRepresentation::SharpJump && !surface) {
        throw std::invalid_argument(
            "sharpJump representation requires a surface constraint.");
    }
    if (config.enforcement == IBMEnforcement::BrinkmanPenalty) {
        if (!body || config.representation != IBMRepresentation::EulerianMask) {
            throw std::invalid_argument(
                "brinkmanPenalty requires body support with eulerianMask.");
        }
        if (!std::isfinite(config.penaltyCoefficient)
            || config.penaltyCoefficient <= 0.0) {
            throw std::invalid_argument(
                "brinkmanPenalty requires a finite penaltyCoefficient > 0.");
        }
    }
    if (config.enforcement == IBMEnforcement::MonolithicKKT && !surface) {
        throw std::invalid_argument(
            "monolithicKKT requires a surface constraint in the current "
            "pressure/velocity/Multiplier formulation.");
    }
    if (config.solidModel == IBMSolidModel::CoupledRigid
        && config.motion != IBMSolidMotion::CoupledRigid) {
        throw std::invalid_argument(
            "solidModel coupledRigid must agree with motion coupledRigid.");
    }
    if (config.motion == IBMSolidMotion::CoupledRigid
        && config.solidModel != IBMSolidModel::CoupledRigid
        && config.solidModel != IBMSolidModel::SelfPropelledRigid) {
        throw std::invalid_argument(
            "motion coupledRigid requires coupledRigid or "
            "selfPropelledRigid solidModel.");
    }
    if (config.solidModel == IBMSolidModel::SelfPropelledRigid
        && config.motion != IBMSolidMotion::CoupledRigid) {
        throw std::invalid_argument(
            "selfPropelledRigid requires motion coupledRigid.");
    }
    validateIBMForcingAlgorithm(config);
}

} // namespace SF::FDM

