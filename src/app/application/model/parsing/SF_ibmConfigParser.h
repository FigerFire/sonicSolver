#pragma once

#include "SF_parserCommon.h"
#include "core/config/SF_configTypes.h"
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace SF::FDM {
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
// 校验已移入 core/config/types/SF_ibmConfigTypes.h（typed 值对象自校验）。
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
} // namespace SF::FDM
