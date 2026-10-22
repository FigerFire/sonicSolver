#pragma once
/// @brief Existing immersed mass/impulse/work relations; no lifecycle or numerical selection.
#include "core/system/SF_formula.h"
#include "core/interfaces/SF_immersedSystem.h"
#include <utility>
namespace SF::IBM::SystemContribution {
namespace Mathematics {
inline Math::Signature signature(std::string_view name,bool surface=true) {
    if(name=="rho")return Math::signature(0,{1,-3,0,0,0,0,0});
    if(name=="rhoU" || name=="rhoU.predictor")return Math::signature(1,{1,-2,-1,0,0,0,0});
    if(name=="rhoE" || name=="rhoE.predictor" || name=="ibm.pCorrection")return Math::signature(0,{1,-1,-2,0,0,0,0});
    for(const auto* id:{"U","U.predictor","U_s","U_s.old","U_b","U_deform"})if(name==id)return Math::signature(1,{0,1,-1,0,0,0,0});
    if(name=="omega_s" || name=="omega_s.old")return Math::signature(1,{0,0,-1,0,0,0,0});
    if(name=="dt")return Math::signature(0,{0,0,1,0,0,0,0});
    if(name=="ibm.forceDensity" || name=="ibm.laggedForceDensity" || name=="lambda_b")return Math::signature(1,{1,-2,-2,0,0,0,0});
    if(name=="Lambda_s")return Math::signature(1,{1,-1,-2,0,0,0,0});
    if(name=="solid.externalForce")return Math::signature(1,{1,1,-2,0,0,0,0});
    if(name=="solid.externalTorque")return Math::signature(1,{1,2,-2,0,0,0,0});
    if(name=="ibm.mask")return Math::signature(0,{});
    if(name=="ibm.penaltyCoefficient" || name=="ibm.gamma")return Math::signature(0,{1,-3,-1,0,0,0,0});
    if(name=="immersed.wallDistance")return Math::signature(0,{0,1,0,0,0,0,0});
    // Domain is a typed product of point coordinates and integration measure;
    // it is never a scalar Field or an operand of ordinary tensor arithmetic.
    if(name=="geometry")return Math::product({Math::signature(1,{0,1,0,0,0,0,0}),Math::signature(0,{0,surface?2.:3.,0,0,0,0,0})});
    throw std::runtime_error("Unknown immersed mathematical symbol: "+std::string(name));
}
inline System::FormulaExpr symbol(const std::string& name,bool surface=true) {return System::FormulaExpr::symbol(name,signature(name,surface));}
inline void check(const std::vector<Math::Signature>& a,std::initializer_list<Math::Signature> types) {
    if(a.size()!=types.size())throw std::runtime_error("Immersed mathematical operator arity mismatch");
    unsigned i=0;for(const auto& t:types)Math::requireCompatible(a[i++],t,"immersed operand");
}
inline System::FormulaExpr op(std::string_view name,std::vector<System::FormulaExpr> args) {
    using S=Math::Signature;using V=std::vector<S>;
    static const Math::OperatorRule rules[]={
        {"interpolate",[](const V& a){check(a,{signature("U")});return a[0];}},
        {"bodyRestrict",[](const V& a){check(a,{signature("U")});return a[0];}},
        {"spread",[](const V& a){check(a,{signature("Lambda_s")});return signature("ibm.forceDensity");}},
        {"bodySpread",[](const V& a){check(a,{signature("lambda_b")});return a[0];}},
        {"midpoint",[](const V& a){if(a.size()!=2)throw std::runtime_error("midpoint needs two values");Math::requireCompatible(a[0],a[1],"midpoint");return a[0];}},
        {"mechanicalWork",[](const V& a){check(a,{signature("rhoU"),signature("U")});return Math::resultSignature(Math::Op::Dot,a);}},
        {"massIncrement",[](const V& a){check(a,{signature("U"),signature("U"),signature("dt")});return signature("ibm.forceDensity");}},
        {"augmentationGradient",[](const V& a){check(a,{signature("U"),signature("U_b"),signature("ibm.gamma")});return Math::resultSignature(Math::Op::Multiply,{a[2],a[0]});}},
        {"translationVelocity",[](const V& a){check(a,{signature("U_s")});return a[0];}},
        {"rotationVelocity",[](const V& a){if(a.size()!=2 || !a[1].tupleSize)throw std::runtime_error("rotation needs omega/domain");Math::requireCompatible(a[0],signature("omega_s"),"angular velocity");return signature("U");}},
        {"rigidVelocity",[](const V& a){check(a,{signature("U_s"),signature("omega_s")});return a[0];}},
        {"diagonalMassResponse",[](const V& a){if(a.size()!=2 || a[1].tupleSize!=2)throw std::runtime_error("mass response needs rho/domain");Math::requireCompatible(a[0],signature("rho"),"response density");const bool surface=(*a[1].tuple)[1].dimensions.exponent[1]==2;Math::requireCompatible(a[1],signature("geometry",surface),"domain measure");return Math::signature(0,{-1,surface?2.:3.,0,0,0,0,0});}},
        {"virtualBodyMass",[](const V& a){if(a.size()!=2 || !a[1].tupleSize)throw std::runtime_error("virtual mass needs density/domain");Math::requireCompatible(a[0],signature("rho"),"virtual mass");return Math::signature(0,{1,0,0,0,0,0,0});}},
        {"virtualBodyMomentum",[](const V& a){if(a.size()!=3 || !a[1].tupleSize)throw std::runtime_error("virtual momentum needs density/domain/velocity");Math::requireCompatible(a[0],signature("rho"),"virtual density");Math::requireCompatible(a[2],signature("U"),"virtual velocity");return Math::signature(1,{1,1,-1,0,0,0,0});}},
        {"virtualBodyInertia",[](const V& a){if(a.size()!=2 || !a[1].tupleSize)throw std::runtime_error("virtual inertia needs density/domain");Math::requireCompatible(a[0],signature("rho"),"virtual density");return Math::signature(2,{1,2,0,0,0,0,0});}},
        {"virtualBodyAngularMomentum",[](const V& a){if(a.size()!=3 || !a[1].tupleSize)throw std::runtime_error("angular momentum needs density/domain/velocity");Math::requireCompatible(a[0],signature("rho"),"virtual density");Math::requireCompatible(a[2],signature("U"),"virtual velocity");return Math::signature(1,{1,2,-1,0,0,0,0});}},
        {"solidMassIncrement",[](const V& a){check(a,{signature("U_s"),signature("U_s.old"),signature("dt")});return signature("solid.externalForce");}},
        {"solidInertiaIncrement",[](const V& a){check(a,{signature("omega_s"),signature("omega_s.old"),signature("dt")});return signature("solid.externalTorque");}},
        {"gyroscopic",[](const V& a){check(a,{signature("omega_s")});return signature("solid.externalTorque");}},
        {"constraintLoad",[](const V& a){check(a,{signature("Lambda_s")});return signature("solid.externalForce");}},
        {"constraintTorque",[](const V& a){check(a,{signature("Lambda_s")});return signature("solid.externalTorque");}}
    };
    for(const auto& r:rules)if(name==r.id)return System::FormulaExpr::op(r,std::move(args));
    return System::FormulaExpr::op(Math::resolveOperator(name),std::move(args));
}
}
}

namespace SF::IBM::SystemContribution {

inline std::vector<System::Equation> mathematics(const FDM::ImmersedAlgorithmDescriptor& selection) {
using namespace SF::System;

using E=FormulaExpr;
const auto s=[&](const std::string& name) { return Mathematics::symbol(name,selection.support==FDM::IBMConstraintSupport::Surface); };
const auto mul=[](E a,E b) {return E::multiply(std::move(a),std::move(b));};
const auto op=[](const std::string& name,std::vector<E> args) {return Mathematics::op(name,std::move(args));};
const bool lagged=selection.enforcement==FDM::IBMEnforcement::ExplicitIBM;
const bool penalty=selection.enforcement==FDM::IBMEnforcement::BrinkmanPenalty;
const bool kkt=selection.enforcement==FDM::IBMEnforcement::MonolithicKKT;
const bool surface=selection.support==FDM::IBMConstraintSupport::Surface;
const bool solid=selection.solid==FDM::IBMSolidModel::SelfPropelledRigid
    || selection.solid==FDM::IBMSolidModel::CoupledRigid;
const std::string lambda=surface ? "Lambda_s" : "lambda_b";
const auto transfer=[&](E value) {return op(surface?"spread":"bodySpread",{std::move(value)});};
const auto rigid= kkt?op("rigidVelocity",{s("U_s"),s("omega_s")}):
    selection.rigidMotionMode==FDM::IBMRigidMotionMode::Rotate
        ?op("rotationVelocity",{s("omega_s"),s("geometry")}):op("translationVelocity",{s("U_s")});
std::vector<SF::System::Equation> equations;
if (lagged) {
    equations.push_back({"ibm.momentum",s("rhoU"),E::add(s("rhoU.predictor"),
        mul(s("dt"),s("ibm.laggedForceDensity"))),{}});
} else if (penalty) {
    equations.push_back({"ibm.momentum",s("rhoU"),E::add(s("rhoU.predictor"),
        mul(s("dt"),op("source",{mul(mul(s("ibm.mask"),s("ibm.penaltyCoefficient")),
            E::subtract(s("U_b"),s("U.predictor")))}))),{}});
} else {
    E increment=op("massIncrement",{s("U"),s("U.predictor"),s("dt")});
    if (kkt) increment=E::add(std::move(increment),op("grad",{s("ibm.pCorrection")}));
    if (selection.algorithm==FDM::IBMForcingAlgorithm::DFMAugmentedLagrangian)
        increment=E::add(std::move(increment),op("augmentationGradient",{s("U"),s("U_b"),s("ibm.gamma")}));
    equations.push_back({"ibm.momentum",E::subtract(std::move(increment),transfer(s(lambda))),
        E::constantValue(0.0),{}});
}
// All current kernels update conservative energy by midpoint mechanical work.
equations.push_back({"ibm.energy",s("rhoE"),E::add(s("rhoE.predictor"),
    op("mechanicalWork",{E::subtract(s("rhoU"),s("rhoU.predictor")),
        op("midpoint",{s("U.predictor"),s("U")})})),{}});
if (!penalty) {
    const auto interpolate=[&](E value) {return op(surface?"interpolate":"bodyRestrict",{std::move(value)});};
    if (lagged) {
        const auto target=solid?E::add(rigid,s("U_deform")):s("U_b");
        equations.push_back({"ibm.multiplier",mul(mul(s("dt"),
            op("diagonalMassResponse",{s("rho"),s("geometry")})),s(lambda)),
            E::subtract(target,interpolate(s("U"))),{}});
    } else {
        E constrained=interpolate(s("U"));
        if (solid) constrained=E::subtract(std::move(constrained),rigid);
        equations.push_back({"ibm.noSlip",std::move(constrained),solid?s("U_deform"):s("U_b"),{}});
    }
}
if (solid && !kkt) {
    // Explicit/fractional self-propulsion uses virtual-fluid generalized mass,
    // not a rigid-body time increment with the configured material mass.
    const bool rotation=selection.rigidMotionMode==FDM::IBMRigidMotionMode::Rotate;
    const auto relative=E::subtract(s(lagged?"U":"U.predictor"),s("U_deform"));
    equations.push_back({"ibm.translation",rotation?s("U_s"):
        mul(op("virtualBodyMass",{s("rho"),s("geometry")}),s("U_s")),
        rotation?s("U_s.old"):E::add(op("virtualBodyMomentum",{s("rho"),s("geometry"),relative}),
            mul(s("dt"),s("solid.externalForce"))),{}});
    equations.push_back({"ibm.rotation",rotation?
        E::op(Math::Op::Dot,{op("virtualBodyInertia",{s("rho"),s("geometry")}),s("omega_s")}):s("omega_s"),
        rotation?E::add(op("virtualBodyAngularMomentum",{s("rho"),s("geometry"),relative}),
            mul(s("dt"),s("solid.externalTorque"))):s("omega_s.old"),{}});
} else if (solid) {
    equations.push_back({"ibm.translation",op("solidMassIncrement",{s("U_s"),s("U_s.old"),s("dt")}),
        E::subtract(s("solid.externalForce"),op("constraintLoad",{s(lambda)})),{}});
    equations.push_back({"ibm.rotation",E::add(op("solidInertiaIncrement",{s("omega_s"),s("omega_s.old"),s("dt")}),
        op("gyroscopic",{s("omega_s")})),E::subtract(s("solid.externalTorque"),
        op("constraintTorque",{s(lambda)})),{}});
}
if (kkt) equations.push_back({"ibm.incompressibility",op("div",{s("U")}),E::constantValue(0.0),{}});
return equations;
}
}
