/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_config.h
/// @brief config 目录总入口：字符串解析和强类型配置校验。
///
/// 强类型枚举和值对象位于 `SF_configTypes.h`。Solver Algorithm 只消费
/// `FDM::SolverConfig`，不直接读取 parser 全局变量。旧 parser 的转换桥接位于
/// `models/initial/SF_legacyConfig.*`。

#include "SF_configTypes.h"
#include "SF_numericsPolicy.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

namespace SF {
namespace FDM {

/// @brief Normalize user-facing config tokens for tolerant parsing.
/// @param value Raw token from `.sf`/TOML-like config or legacy caller.
/// @return Lowercase token with quotes, underscores, dashes, and whitespace removed.
inline std::string normalizeToken(std::string value) {
    value.erase(std::remove(value.begin(), value.end(), '"'), value.end());
    value.erase(std::remove_if(value.begin(), value.end(),
                               [](unsigned char c) {
                                   return c == '_' || c == '-' || std::isspace(c);
                               }),
                value.end());
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return (char)std::tolower(c); });
    return value;
}

/// @brief Split a source-list style string.
/// @param value Text using comma, semicolon, plus, or whitespace separators.
/// @return Ordered non-empty tokens; tokens are not normalized here.
inline std::vector<std::string> splitList(const std::string& value) {
    std::vector<std::string> out;
    std::string token;
    for (char c : value) {
        if (c == ',' || c == ';' || c == '+' || std::isspace((unsigned char)c)) {
            if (!token.empty()) {
                out.push_back(token);
                token.clear();
            }
        } else {
            token.push_back(c);
        }
    }
    if (!token.empty()) out.push_back(token);
    return out;
}

/// @brief Parse a convective scheme name.
/// @param value User-facing string such as `"WENO3"`, `"WENO5"`, `"TENO5"` or `"WENO7"`.
/// @return Strongly typed scheme; empty input keeps `WENO5`.
inline ConvectionScheme parseConvectionScheme(const std::string& value) {
    std::string t = normalizeToken(value);
    if (t.empty()) return ConvectionScheme::WENO5;
    if (t == "weno3") return ConvectionScheme::WENO3;
    if (t == "weno5") return ConvectionScheme::WENO5;
    if (t == "teno5") return ConvectionScheme::TENO5;
    if (t == "weno7") return ConvectionScheme::WENO7;
    throw std::invalid_argument(
        "Unsupported convection scheme '" + value
        + "'. Supported schemes: WENO3, WENO5, TENO5, WENO7. "
          "Recommended pairings: WENO5/TENO5 with StegerWarming or "
          "LaxFriedrichs for smooth shock-capturing cases; WENO3 with "
          "LaxFriedrichs for IBM/ILW debug runs; WENO7 only when ghost depth "
          "and boundary closure are known to be sufficient.");
}

/// @brief Parse the governing-equation formulation.
/// @param value User-facing string such as `"conservativeFluxDifference"`.
/// @return Strongly typed formulation; empty input keeps the conservative default.
inline EquationFormulation parseEquationFormulation(const std::string& value) {
    std::string t = normalizeToken(value);
    if (t.empty() || t == "conservative" || t == "conservativefd"
        || t == "conservativefinitedifference"
        || t == "conservativefluxdifference"
        || t == "fluxdifference" || t == "fluxdifferencing") {
        return EquationFormulation::ConservativeFluxDifference;
    }
    if (t == "primitive" || t == "primitivedifferential"
        || t == "nonconservative" || t == "differential") {
        return EquationFormulation::PrimitiveDifferential;
    }
    throw std::invalid_argument(
        "Unsupported numerics formulation '" + value
        + "'. Use conservativeFluxDifference for density-based compressible flow.");
}

/// @brief Parse the variable family reconstructed at faces.
/// @param value User-facing string such as `"characteristic"`.
/// @return Strongly typed reconstruction variable; empty input keeps characteristic.
inline ReconstructionVariable parseReconstructionVariable(const std::string& value) {
    std::string t = normalizeToken(value);
    if (t.empty() || t == "characteristic" || t == "char") {
        return ReconstructionVariable::Characteristic;
    }
    if (t == "conservative" || t == "u") return ReconstructionVariable::Conservative;
    if (t == "primitive" || t == "prim") return ReconstructionVariable::Primitive;
    throw std::invalid_argument(
        "Unsupported reconstruction variable '" + value
        + "'. Current WENO flux assembly supports characteristic reconstruction.");
}

/// @brief Parse the common-face interface flux policy.
/// @param value User-facing string such as `"shared"` or `"sharedInterfaceFlux"`.
/// @return Strongly typed interface policy.
inline InterfaceFluxPolicy parseInterfaceFluxPolicy(const std::string& value) {
    std::string t = normalizeToken(value);
    if (t.empty() || t == "shared" || t == "sharedflux"
        || t == "sharedinterfaceflux" || t == "commonfaceflux") {
        return InterfaceFluxPolicy::SharedInterfaceFlux;
    }
    throw std::invalid_argument(
        "Unsupported interface flux policy '" + value
        + "'. Current conservative solver requires sharedInterfaceFlux.");
}

/// @brief Parse an inviscid flux method name.
/// @param value User-facing string such as `"StegerWarming"`, `"Rusanov"` or `"Roe"`.
/// @return Strongly typed flux method; empty input keeps `StegerWarming`.
inline FluxSplitter parseFluxSplitter(const std::string& value) {
    std::string t = normalizeToken(value);
    if (t.empty()) return FluxSplitter::StegerWarming;
    if (t == "rusanov" || t == "localaxfriedrichs" || t == "llf") return FluxSplitter::Rusanov;
    if (t == "laxfriedrichs" || t == "lf" || t == "lxf" || t == "laxf") return FluxSplitter::LaxFriedrichs;
    if (t == "roe") return FluxSplitter::Roe;
    if (t == "laxwendroff" || t == "lw") return FluxSplitter::LaxWendroff;
    if (t == "stegerwarming" || t == "sw") return FluxSplitter::StegerWarming;
    throw std::invalid_argument(
        "Unsupported inviscid flux method '" + value
        + "'. Supported methods: StegerWarming, Rusanov, LaxFriedrichs, Roe, "
          "LaxWendroff. Recommended pairings: WENO5/TENO5+StegerWarming for "
          "regular cases, WENO3+LaxFriedrichs for IBM/ILW debug, Roe for "
          "direct face-flux checks.");
}

/// @brief Parse the IBM near-wall WENO closure mode.
/// @param value User-facing string such as `"Downgrade"`, `"LowOrder"` or `"ILW"`.
/// @return Strongly typed IBM boundary scheme; empty input keeps low-order.
inline IBMBoundaryScheme parseIBMBoundaryScheme(const std::string& value) {
    std::string t = normalizeToken(value);
    if (t.empty() || t == "downgrade" || t == "loworder"
        || t == "firstorder") {
        return IBMBoundaryScheme::LowOrder;
    }
    if (t == "ilw" || t == "highorder" || t == "highorderilw") {
        return IBMBoundaryScheme::ILW;
    }
    throw std::invalid_argument(
        "Unsupported IBM WENO closure '" + value
        + "'. Supported closures: Downgrade/LowOrder or ILW. "
          "Use ILW only when [numerics].ILW is configured with the requested "
          "WENO order; use Downgrade for explicit low-order IBM debugging.");
}

/// @brief 解析 IBM 主方法；未实现的 direct/continuous forcing 不会被替换。
inline IBMMethod parseIBMMethod(const std::string& value) {
    const std::string token = normalizeToken(value);
    if (token.empty() || token == "ghost" || token == "ghostcell") {
        return IBMMethod::Ghost;
    }
    if (token == "variationalforcing" || token == "variationaldlm"
        || token == "constraintforcing") {
        return IBMMethod::VariationalForcing;
    }
    throw std::invalid_argument(
        "Unsupported IBM type '" + value
        + "'. Implemented types: ghost, variationalForcing. "
          "directForcing/continuousForcing are not aliases and are not "
          "silently substituted.");
}

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

/// @brief 解析论文 Algorithms 1--7 对应的显式算法名。
inline IBMForcingAlgorithm parseIBMForcingAlgorithm(
        const std::string& value) {
    const std::string token = normalizeToken(value);
    if (token == "peskin" || token == "peskinoriginal"
        || token == "originalibm") {
        return IBMForcingAlgorithm::PeskinOriginal;
    }
    if (token == "dfmexplicitselfpropelled"
        || token == "explicitselfpropelled") {
        return IBMForcingAlgorithm::DFMExplicitSelfPropelled;
    }
    if (token == "dfmfractionalstepselfpropelled"
        || token == "dfmftsselfpropelled"
        || token == "fractionalselfpropelled") {
        return IBMForcingAlgorithm::DFMFractionalStepSelfPropelled;
    }
    if (token == "dfmfractionalstepprescribed"
        || token == "dfmftsprescribed"
        || token == "fractionaldlm") {
        return IBMForcingAlgorithm::DFMFractionalStepPrescribed;
    }
    if (token == "dfmimplicitprescribed"
        || token == "fullyimplicitprescribed") {
        return IBMForcingAlgorithm::DFMImplicitPrescribed;
    }
    if (token == "dfmimplicitselfpropelled"
        || token == "fullyimplicitselfpropelled") {
        return IBMForcingAlgorithm::DFMImplicitSelfPropelled;
    }
    if (token == "dfmaugmentedlagrangian"
        || token == "augmentedlagrangian"
        || token == "augmentedkkt") {
        return IBMForcingAlgorithm::DFMAugmentedLagrangian;
    }
    if (token == "velocityforcingfts" || token == "directforcingfts") {
        return IBMForcingAlgorithm::VelocityForcingFTS;
    }
    if (token == "velocityforcingbp" || token == "brinkman"
        || token == "brinkmanpenalty") {
        return IBMForcingAlgorithm::VelocityForcingBP;
    }
    throw std::invalid_argument(
        "IBM forcing algorithm must be peskinOriginal, "
        "dfmExplicitSelfPropelled, dfmFractionalStepSelfPropelled, "
        "dfmFractionalStepPrescribed, dfmImplicitPrescribed, "
        "dfmImplicitSelfPropelled, dfmAugmentedLagrangian, "
        "velocityForcingFTS, or "
        "velocityForcingBP; got '" + value + "'.");
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

inline IBMConstraintDomain parseIBMConstraintDomain(
        const std::string& value) {
    const std::string token = normalizeToken(value);
    if (token == "volume" || token == "bodyvolume") {
        return IBMConstraintDomain::Volume;
    }
    if (token == "surface" || token == "boundarysurface") {
        return IBMConstraintDomain::Surface;
    }
    throw std::invalid_argument(
        "IBM forcing constraintDomain must be volume or surface; got '"
        + value + "'.");
}

/// @brief 解析正交的浸没约束支持域。
inline IBMConstraintSupport parseIBMConstraintSupport(
        const std::string& value) {
    const std::string token = normalizeToken(value);
    if (token == "body" || token == "volume" || token == "bodyvolume") {
        return IBMConstraintSupport::Body;
    }
    if (token == "surface" || token == "boundarysurface") {
        return IBMConstraintSupport::Surface;
    }
    if (token == "surfaceandbody" || token == "bodysurface"
        || token == "both" || token == "volumeandsurface") {
        return IBMConstraintSupport::SurfaceAndBody;
    }
    throw std::invalid_argument(
        "IBM forcing constraintSupport must be body, surface, or "
        "surfaceAndBody; got '" + value + "'.");
}

/// @brief 解析浸没力的空间表示。
inline IBMRepresentation parseIBMRepresentation(
        const std::string& value) {
    const std::string token = normalizeToken(value);
    if (token == "diffusekernel" || token == "kernel"
        || token == "lagrangiankernel") {
        return IBMRepresentation::DiffuseKernel;
    }
    if (token == "eulerianmask" || token == "eulerianbody"
        || token == "volumemask" || token == "bodymask") {
        return IBMRepresentation::EulerianMask;
    }
    if (token == "sharpjump" || token == "sharpinterface"
        || token == "jump") {
        return IBMRepresentation::SharpJump;
    }
    throw std::invalid_argument(
        "IBM forcing representation must be diffuseKernel, eulerianMask, "
        "or sharpJump; got '" + value + "'.");
}

/// @brief 解析浸没约束的代数/时间强制策略。
inline IBMEnforcement parseIBMEnforcement(const std::string& value) {
    const std::string token = normalizeToken(value);
    if (token == "explicitibm" || token == "classicalibm"
        || token == "laggedibm" || token == "explicitmultiplier") {
        return IBMEnforcement::ExplicitIBM;
    }
    if (token == "fractionaldlm" || token == "fictitiousdomain"
        || token == "fractionalprojection" || token == "dlm") {
        return IBMEnforcement::FractionalDLM;
    }
    if (token == "velocityforcing" || token == "directforcing"
        || token == "fts") {
        return IBMEnforcement::VelocityForcing;
    }
    if (token == "brinkman" || token == "brinkmanpenalty"
        || token == "penalty") {
        return IBMEnforcement::BrinkmanPenalty;
    }
    if (token == "monolithickkt" || token == "fullyimplicit"
        || token == "implicitkkt" || token == "kkt") {
        return IBMEnforcement::MonolithicKKT;
    }
    throw std::invalid_argument(
        "IBM forcing enforcement must be explicitIBM, fractionalDLM, "
        "velocityForcing, brinkmanPenalty, or monolithicKKT; got '"
        + value + "'.");
}

/// @brief 解析固体未知量模型。
inline IBMSolidModel parseIBMSolidModel(const std::string& value) {
    const std::string token = normalizeToken(value);
    if (token == "prescribed" || token == "prescribedrigid") {
        return IBMSolidModel::Prescribed;
    }
    if (token == "coupledrigid" || token == "rigid6dof"
        || token == "rigid") {
        return IBMSolidModel::CoupledRigid;
    }
    if (token == "selfpropelled" || token == "selfpropelledrigid") {
        return IBMSolidModel::SelfPropelledRigid;
    }
    if (token == "deformable" || token == "elastic") {
        return IBMSolidModel::Deformable;
    }
    throw std::invalid_argument(
        "IBM forcing solidModel must be prescribed, coupledRigid, "
        "selfPropelledRigid, or deformable; got '" + value + "'.");
}

/// @brief 校验 algorithm 与空间/时间离散轴是否描述同一方法。
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

inline IBMConstraintCoupling parseIBMConstraintCoupling(
        const std::string& value) {
    const std::string token = normalizeToken(value);
    if (token == "incrementalprojection") {
        return IBMConstraintCoupling::IncrementalProjection;
    }
    if (token == "monolithickkt" || token == "kkt") {
        return IBMConstraintCoupling::MonolithicKKT;
    }
    throw std::invalid_argument(
        "IBM forcing coupling must be incrementalProjection or "
        "monolithicKKT; got '" + value + "'.");
}

inline IBMSolidMotion parseIBMSolidMotion(const std::string& value) {
    const std::string token = normalizeToken(value);
    if (token == "prescribedrigid") {
        return IBMSolidMotion::PrescribedRigid;
    }
    if (token == "coupledrigid") {
        return IBMSolidMotion::CoupledRigid;
    }
    throw std::invalid_argument(
        "IBM forcing motion must be prescribedRigid or coupledRigid; got '"
        + value + "'.");
}

/// @brief 解析刚体运动自由度；motivation 是平动模式的用户关键字。
inline IBMRigidMotionMode parseIBMRigidMotionMode(
        const std::string& value) {
    const std::string token = normalizeToken(value);
    if (token == "motivation" || token == "translation"
        || token == "translate") {
        return IBMRigidMotionMode::Motivation;
    }
    if (token == "rotate" || token == "rotation") {
        return IBMRigidMotionMode::Rotate;
    }
    throw std::invalid_argument(
        "IBM forcing motionMode must be motivation or rotate; got '"
        + value + "'.");
}

inline IBMEnergyCoupling parseIBMEnergyCoupling(
        const std::string& value) {
    const std::string token = normalizeToken(value);
    if (token == "mechanicalwork") {
        return IBMEnergyCoupling::MechanicalWork;
    }
    throw std::invalid_argument(
        "IBM forcing energyCoupling must be mechanicalWork; got '"
        + value + "'.");
}

/// @brief Parse a time-integration scheme name.
/// @param value User-facing string such as `"Euler"`, `"SSP_RK3"` or `"RK4"`.
/// @return Strongly typed time scheme; empty input keeps `RK4`.
inline TimeScheme parseTimeScheme(const std::string& value) {
    std::string t = normalizeToken(value);
    if (t.empty()) return TimeScheme::RK4;
    if (t == "euler") return TimeScheme::Euler;
    if (t == "ssprk3" || t == "rk3") return TimeScheme::SSPRK3;
    if (t == "rk4" || t == "rungekutta4" || t == "classicrk4") {
        return TimeScheme::RK4;
    }
    throw std::invalid_argument(
        "Unsupported time scheme '" + value
        + "'. Supported schemes: Euler, SSP_RK3, RK4.");
}

/// @brief Parse a viscous central-difference scheme name.
/// @param value User-facing string such as `"CENTRAL2"` or `"CENTRAL4"`.
/// @return Strongly typed viscous scheme; empty input keeps `CENTRAL2`.
inline ViscousScheme parseViscousScheme(const std::string& value) {
    std::string t = normalizeToken(value);
    if (t.empty() || t == "central2" || t == "cd2") return ViscousScheme::Central2;
    if (t == "central4" || t == "cd4") return ViscousScheme::Central4;
    throw std::invalid_argument(
        "Unsupported viscous scheme '" + value
        + "'. Supported schemes: CENTRAL2, CENTRAL4. Recommended pairing: "
          "CENTRAL2 for current viscous/turbulence validation; CENTRAL4 only "
          "after boundary stencils are verified.");
}

/// @brief Parse source-term family names.
/// @param value Source list such as `"Gravity+MRF"` or `"None"`.
/// @return Ordered source kinds; disabled/none tokens are omitted.
inline std::vector<SourceKind> parseSourceKinds(const std::string& value) {
    std::vector<SourceKind> kinds;
    std::string sourceText = value.empty() ? "None" : value;
    for (const std::string& raw : splitList(sourceText)) {
        std::string t = normalizeToken(raw);
        if (t.empty() || t == "none" || t == "off" || t == "false") continue;
        if (t == "gravity") kinds.push_back(SourceKind::Gravity);
        else if (t == "mrf" || t == "rotating" || t == "rotation") kinds.push_back(SourceKind::MRF);
        else if (t == "wallheat" || t == "wallheatsource"
                 || t == "wallheatflux" || t == "heat" || t == "heater") {
            kinds.push_back(SourceKind::WallHeat);
        }
        else {
            throw std::invalid_argument(
                "Unsupported source term '" + raw
                + "'. Supported source terms: None, Gravity, MRF, WallHeat.");
        }
    }
    return kinds;
}

/// @brief 解析湍流大类字符串。
/// @param value 用户配置中的 type，例如 `"RAS"`、`"LES"`、`"DNS"`。
/// @return 强类型湍流大类；空值返回 `None`。
inline TurbulenceFamily parseTurbulenceFamily(const std::string& value) {
    std::string t = normalizeToken(value);
    if (t.empty() || t == "none" || t == "off" || t == "false"
        || t == "laminar") {
        return TurbulenceFamily::None;
    }
    if (t == "ras" || t == "rans") return TurbulenceFamily::RAS;
    if (t == "les") return TurbulenceFamily::LES;
    if (t == "dns") return TurbulenceFamily::DNS;
    throw std::invalid_argument(
        "Unsupported turbulence family '" + value
        + "'. Supported families: None, RAS, LES, DNS. Recommended "
          "pairings: RAS+kEpsilon/kOmegaSST, LES+Smagorinsky, DNS+DNS.");
}

/// @brief 解析具体湍流模型字符串。
/// @param value 用户配置中的 value/model，例如 `"kEpsilon"`、`"Smagorinsky"`。
/// @return 强类型模型标识；空值返回 `None`。
inline TurbulenceModelKind parseTurbulenceModelKind(const std::string& value) {
    std::string t = normalizeToken(value);
    if (t.empty() || t == "none" || t == "off" || t == "false"
        || t == "laminar") {
        return TurbulenceModelKind::None;
    }
    if (t == "kepsilon" || t == "kep" || t == "ke") return TurbulenceModelKind::kEpsilon;
    if (t == "komegasst" || t == "sst") return TurbulenceModelKind::kOmegaSST;
    if (t == "smagorinsky" || t == "smago") return TurbulenceModelKind::Smagorinsky;
    if (t == "dns") return TurbulenceModelKind::DNS;
    throw std::invalid_argument(
        "Unsupported turbulence model '" + value
        + "'. Supported models: None, kEpsilon, kOmegaSST, Smagorinsky, DNS. "
          "Recommended pairings: RAS+kEpsilon/kOmegaSST, LES+Smagorinsky, "
          "DNS+DNS.");
}

/// @brief Convert a convective scheme to the canonical config/log token.
/// @param scheme Strongly typed scheme.
/// @return Stable string used in logs and dispatch bridges.
inline const char* toString(ConvectionScheme scheme) {
    switch (scheme) {
        case ConvectionScheme::WENO3: return "WENO3";
        case ConvectionScheme::WENO5: return "WENO5";
        case ConvectionScheme::TENO5: return "TENO5";
        case ConvectionScheme::WENO7: return "WENO7";
    }
    return "WENO5";
}

/// @brief Convert formulation to the canonical config/log token.
inline const char* toString(EquationFormulation formulation) {
    switch (formulation) {
        case EquationFormulation::ConservativeFluxDifference:
            return "conservativeFluxDifference";
        case EquationFormulation::PrimitiveDifferential:
            return "primitiveDifferential";
    }
    return "conservativeFluxDifference";
}

/// @brief Convert reconstruction variable to the canonical config/log token.
inline const char* toString(ReconstructionVariable variable) {
    switch (variable) {
        case ReconstructionVariable::Characteristic: return "characteristic";
        case ReconstructionVariable::Conservative: return "conservative";
        case ReconstructionVariable::Primitive: return "primitive";
    }
    return "characteristic";
}

/// @brief Convert interface flux policy to the canonical config/log token.
inline const char* toString(InterfaceFluxPolicy policy) {
    switch (policy) {
        case InterfaceFluxPolicy::SharedInterfaceFlux:
            return "sharedInterfaceFlux";
    }
    return "sharedInterfaceFlux";
}

/// @brief Convert a flux method to the canonical config/log token.
/// @param splitter Strongly typed flux method.
/// @return Stable string used in logs and dispatch bridges.
inline const char* toString(FluxSplitter splitter) {
    switch (splitter) {
        case FluxSplitter::StegerWarming: return "StegerWarming";
        case FluxSplitter::Rusanov: return "Rusanov";
        case FluxSplitter::LaxFriedrichs: return "LaxFriedrichs";
        case FluxSplitter::Roe: return "Roe";
        case FluxSplitter::LaxWendroff: return "LaxWendroff";
    }
    return "unknown";
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

/// @brief Convert a time scheme to the canonical config/log token.
/// @param scheme Strongly typed time scheme.
/// @return Stable string used in logs and dispatch bridges.
inline const char* toString(TimeScheme scheme) {
    switch (scheme) {
        case TimeScheme::Euler: return "Euler";
        case TimeScheme::SSPRK3: return "SSP_RK3";
        case TimeScheme::RK4: return "RK4";
    }
    return "RK4";
}

/// @brief Return the number of explicit stages owned by a time scheme.
/// @throws std::invalid_argument if the enum value is not implemented.
inline int explicitStageCount(TimeScheme scheme) {
    switch (scheme) {
        case TimeScheme::Euler: return 1;
        case TimeScheme::SSPRK3: return 3;
        case TimeScheme::RK4: return 4;
    }
    throw std::invalid_argument("Unsupported explicit time scheme.");
}

/// @brief Convert a viscous scheme to the canonical config/log token.
/// @param scheme Strongly typed viscous scheme.
/// @return Stable string used in logs and dispatch bridges.
inline const char* toString(ViscousScheme scheme) {
    switch (scheme) {
        case ViscousScheme::Central2: return "CENTRAL2";
        case ViscousScheme::Central4: return "CENTRAL4";
    }
    return "CENTRAL2";
}

/// @brief 转换湍流大类到稳定日志字符串。
/// @param family 强类型湍流大类。
/// @return 稳定字符串。
inline const char* toString(TurbulenceFamily family) {
    switch (family) {
        case TurbulenceFamily::None: return "None";
        case TurbulenceFamily::RAS: return "RAS";
        case TurbulenceFamily::LES: return "LES";
        case TurbulenceFamily::DNS: return "DNS";
    }
    return "None";
}

/// @brief 转换具体湍流模型到稳定日志字符串。
/// @param model 强类型模型枚举。
/// @return 稳定字符串。
inline const char* toString(TurbulenceModelKind model) {
    switch (model) {
        case TurbulenceModelKind::None: return "None";
        case TurbulenceModelKind::kEpsilon: return "kEpsilon";
        case TurbulenceModelKind::kOmegaSST: return "kOmegaSST";
        case TurbulenceModelKind::Smagorinsky: return "Smagorinsky";
        case TurbulenceModelKind::DNS: return "DNS";
    }
    return "None";
}

/// @brief Convert source kinds to a compact source-list token.
/// @param kinds Ordered source kinds.
/// @return `"None"` when empty, otherwise names joined with `+`.
inline std::string toString(const std::vector<SourceKind>& kinds) {
    if (kinds.empty()) return "None";

    std::ostringstream out;
    for (size_t i = 0; i < kinds.size(); ++i) {
        if (i > 0) out << "+";
        switch (kinds[i]) {
            case SourceKind::Gravity: out << "Gravity"; break;
            case SourceKind::MRF: out << "MRF"; break;
            case SourceKind::WallHeat: out << "WallHeat"; break;
        }
    }
    return out.str();
}

/// @brief 校验线性求解配置，不允许无效容差或隐式算法替换。
inline void validateLinearSolverConfig(const LinearSolverConfig& config,
                                       const std::string& context) {
    if (config.maxIterations <= 0 || config.krylovDimension <= 0
        || !std::isfinite(config.relativeTolerance)
        || config.relativeTolerance <= 0.0
        || !std::isfinite(config.absoluteTolerance)
        || config.absoluteTolerance < 0.0
        || config.structureRebuildInterval < 0
        || config.preconditionerRefreshInterval < 0
        || config.iluLevelOfFill < 0
        || config.amgSweeps <= 0 || config.amgMaxLevels <= 0
        || config.amgCoarsenType < 0 || config.amgRelaxType < 0
        || !std::isfinite(config.amgStrongThreshold)
        || config.amgStrongThreshold <= 0 || config.amgStrongThreshold > 1
        || config.schurDenseLimit <= 0) {
        throw std::invalid_argument(
            context + ": invalid Krylov iteration/tolerance controls.");
    }
    if(config.method==KrylovMethod::PCG
       &&config.preconditioner!=LinearPreconditioner::BoomerAMG) {
        throw std::invalid_argument(
            context+": PCG requires the symmetric boomerAMG preconditioner.");
    }
    if(config.method==KrylovMethod::PCG
       &&config.equilibration!=LinearEquilibration::None)
        throw std::invalid_argument(context+": rowMax scaling is not symmetric; use flexGMRES.");
    if(config.preconditioner==LinearPreconditioner::ILU
       &&config.iluType!=0 && config.iluType!=1
       &&config.iluType!=10 && config.iluType!=11
       &&config.iluType!=20 && config.iluType!=21
       &&config.iluType!=30 && config.iluType!=31
       &&config.iluType!=40 && config.iluType!=41
       &&config.iluType!=50) {
        throw std::invalid_argument(
            context+": unsupported HYPRE ILU type.");
    }
}

/// @brief 校验 SIMPLE/PISO/PIMPLE 与压力参考配置。
inline void validateSolverProperties(const SolverPropertiesConfig& config) {
    if (config.outerCorrectors <= 0 || config.pressureCorrectors <= 0
        || config.nonOrthogonalCorrectors < 0
        || !std::isfinite(config.momentumRelaxation)
        || config.momentumRelaxation <= 0.0
        || config.momentumRelaxation > 1.0
        || !std::isfinite(config.pressureRelaxation)
        || config.pressureRelaxation <= 0.0
        || config.pressureRelaxation > 1.0
        || !std::isfinite(config.phaseSourceCfl)
        || config.phaseSourceCfl <= 0.0
        || config.phaseSourceCfl > 1.0) {
        throw std::invalid_argument(
            "solverProperties has invalid corrector or relaxation controls.");
    }
    if (config.referenceCell < 0
        || !std::isfinite(config.referencePressure)
        || config.referencePressure <= 0.0) {
        throw std::invalid_argument(
            "pressure workflow requires explicit non-negative referenceCell "
            "and finite referencePressure > 0.");
    }
    validateLinearSolverConfig(config.pressure, "pressure linear solver");
    validateLinearSolverConfig(config.momentum, "momentum linear solver");
    validateLinearSolverConfig(config.energy, "energy linear solver");
    validateLinearSolverConfig(config.turbulence, "turbulence linear solver");
}

/// @brief Fail fast for declared numerical modes that do not have an implementation.
/// @param numerics Parsed numerical-method configuration.
inline void validateNumericsConfig(const NumericsConfig& numerics) {
    if (!std::isfinite(numerics.cfl) || numerics.cfl <= 0.0) {
        throw std::invalid_argument(
            "CFL/maxCo must be a finite positive value.");
    }
    if (!std::isfinite(numerics.maxDeltaT)
        || numerics.maxDeltaT <= 0.0) {
        throw std::invalid_argument(
            "maxDeltaT must be a finite positive value.");
    }
    if (!std::isfinite(numerics.idealGasGamma)
        || numerics.idealGasGamma <= 1.0
        || !std::isfinite(numerics.idealGasConstant)
        || numerics.idealGasConstant <= 0.0) {
        throw std::invalid_argument(
            "legacy ideal-gas numerics require gamma > 1 and R > 0.");
    }
    if (numerics.formulation
        != EquationFormulation::ConservativeFluxDifference) {
        throw std::invalid_argument(
            "primitive/non-conservative differential formulation is not "
            "supported for the density-based compressible solver; use "
            "conservativeFluxDifference.");
    }
    if (numerics.reconstruction
        != ReconstructionVariable::Characteristic) {
        throw std::invalid_argument(
            "requested face reconstruction '"
            + std::string(toString(numerics.reconstruction))
            + "' is not implemented yet; current WENO flux assembly supports "
            "characteristic reconstruction.");
    }
}

} // namespace FDM
} // namespace SF
