#pragma once
/// @brief Existing immersed mass/impulse/work relations; no lifecycle or numerical selection.
#include "core/system/SF_formula.h"
#include "core/interfaces/SF_immersedSystem.h"
#include <utility>
namespace SF::IBM::SystemContribution {
inline std::vector<System::Equation> mathematics(const FDM::ImmersedAlgorithmDescriptor& selection) {
using namespace SF::System;

using E=FormulaExpr;
const auto s=[](const std::string& name) { return E::symbol(name); };
const auto mul=[](E a,E b) {return E::multiply(std::move(a),std::move(b));};
const auto op=[](const std::string& name,std::vector<E> args) {return E::op(name,std::move(args));};
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
        mul(s("dt"),s("ibm.laggedForceDensity"))),{},true});
} else if (penalty) {
    equations.push_back({"ibm.momentum",s("rhoU"),E::add(s("rhoU.predictor"),
        mul(s("dt"),op("source",{mul(mul(s("ibm.mask"),s("ibm.penaltyCoefficient")),
            E::subtract(s("U_b"),s("U.predictor")))}))),{},true});
} else {
    E increment=op("massIncrement",{s("U"),s("U.predictor"),s("dt")});
    if (kkt) increment=E::add(std::move(increment),op("grad",{s("ibm.pCorrection")}));
    if (selection.algorithm==FDM::IBMForcingAlgorithm::DFMAugmentedLagrangian)
        increment=E::add(std::move(increment),op("augmentationGradient",{s("U"),s("U_b"),s("ibm.gamma")}));
    equations.push_back({"ibm.momentum",E::subtract(std::move(increment),transfer(s(lambda))),
        E::constantValue(0.0),{},true});
}
// All current kernels update conservative energy by midpoint mechanical work.
equations.push_back({"ibm.energy",s("rhoE"),E::add(s("rhoE.predictor"),
    op("mechanicalWork",{E::subtract(s("rhoU"),s("rhoU.predictor")),
        op("midpoint",{s("U.predictor"),s("U")})})),{},true});
if (!penalty) {
    const auto interpolate=[&](E value) {return op(surface?"interpolate":"bodyRestrict",{std::move(value)});};
    if (lagged) {
        const auto target=solid?E::add(rigid,s("U_deform")):s("U_b");
        equations.push_back({"ibm.multiplier",mul(mul(s("dt"),
            op("diagonalMassResponse",{s("rho"),s("geometry")})),s(lambda)),
            E::subtract(target,interpolate(s("U"))),{},true});
    } else {
        E constrained=interpolate(s("U"));
        if (solid) constrained=E::subtract(std::move(constrained),rigid);
        equations.push_back({"ibm.noSlip",std::move(constrained),solid?s("U_deform"):s("U_b"),{},true});
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
            mul(s("dt"),s("solid.externalForce"))),{},true});
    equations.push_back({"ibm.rotation",rotation?
        mul(op("virtualBodyInertia",{s("rho"),s("geometry")}),s("omega_s")):s("omega_s"),
        rotation?E::add(op("virtualBodyAngularMomentum",{s("rho"),s("geometry"),relative}),
            mul(s("dt"),s("solid.externalTorque"))):s("omega_s.old"),{},true});
} else if (solid) {
    equations.push_back({"ibm.translation",op("solidMassIncrement",{s("U_s"),s("U_s.old"),s("dt")}),
        E::subtract(s("solid.externalForce"),op("constraintLoad",{s(lambda)})),{},true});
    equations.push_back({"ibm.rotation",E::add(op("solidInertiaIncrement",{s("omega_s"),s("omega_s.old"),s("dt")}),
        op("gyroscopic",{s("omega_s")})),E::subtract(s("solid.externalTorque"),
        op("constraintTorque",{s(lambda)})),{},true});
}
if (kkt) equations.push_back({"ibm.incompressibility",op("div",{s("U")}),E::constantValue(0.0),{},true});
return equations;
}
}
