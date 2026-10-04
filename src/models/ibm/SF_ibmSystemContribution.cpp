/// @file SF_ibmSystemContribution.cpp
/// @brief IBM 模型注册 boundary 或原生 impulse/projection/block 数学，不拥有 lifecycle。
#include "SF_ibmSystemContribution.h"
#include "SF_immersedMathematics.h"
#include "SF_immersedSystem.h"
#include <utility>

namespace SF::IBM::SystemContribution {
using namespace SF::System;
void contribute(System::SystemContribution& system,
                const FDM::ImmersedAlgorithmDescriptor& immersed) {
    system.recordContribution("model.ibm."+immersed.id,"immersed-boundary contribution");
    if (immersed.enforcement==FDM::IBMEnforcement::GhostCell) {
        system.requireProvider("ibm.boundary","bind stage-time ghost/ILW stencil closure");
        system.addBoundaryClosure({"immersed.ghost","ibm.boundary",{"conservative","geometry","classification"},
            {"conservative.ghost"},true,
            {"physical boundary","halo COPY","ghost/ILW reconstruction","ghost publication","halo COPY","spatial operator"},{}});
        return;
    }
    system.requireProvider("ibm.constraint","execute compiled immersed impulse/projection/block");
    const auto state=[&](std::string id,int components,VariableLocation location,
                         OwnershipKind owner,StorageBinding binding,std::string storage,StateRole role) {
        StateSymbol value;
        value.id=id;value.name=id;value.components=components;
        value.shape=components==1?ValueShape::Scalar:ValueShape::Vector;
        value.location=location;value.ownership=owner;value.storageBinding=binding;
        value.storageKey=std::move(storage);value.role=role;value.nameSpace="immersed";
        value.initializationRequired=false;value.boundaryRequired=false;
        value.restartEligible=false;value.outputEligible=false;
        system.addState(std::move(value));
    };
    state("ibm.forceDensity",3,VariableLocation::EulerianCell,OwnershipKind::EulerianGlobalDof,
        StorageBinding::SpecializedExecutor,"ibm.forceDensity",StateRole::Derived);
    if (immersed.enforcement==FDM::IBMEnforcement::ExplicitIBM)
        state("ibm.laggedForceDensity",3,VariableLocation::EulerianCell,OwnershipKind::EulerianGlobalDof,
            StorageBinding::SpecializedExecutor,"ibm.laggedForceDensity",StateRole::Algebraic);
    const bool surface=immersed.support==FDM::IBMConstraintSupport::Surface;
    const bool multiplier=immersed.enforcement!=FDM::IBMEnforcement::BrinkmanPenalty;
    const std::string lambda=surface?"Lambda_s":"lambda_b";
    if (multiplier) {
        // Surface lambda is local solver workspace; existing kernels do not retain
        // a persistent surface array. The body multiplier aliases the force array.
        state(lambda,3,surface?VariableLocation::SurfaceConstraint:VariableLocation::BodyConstraint,
            OwnershipKind::ConstraintGlobalDof,
            surface?StorageBinding::TransientWorkspace:StorageBinding::SpecializedExecutor,
            surface?"ibm.surfaceMultiplier":"ibm.forceDensity",StateRole::Multiplier);
    }
    const bool coupled=immersed.solid==FDM::IBMSolidModel::SelfPropelledRigid
        || immersed.solid==FDM::IBMSolidModel::CoupledRigid;
    if (coupled) for (const auto* id:{"U_s","omega_s"})
        state(id,3,VariableLocation::SolidGlobal,OwnershipKind::SolidGlobalDof,
            StorageBinding::SpecializedExecutor,std::string("ibm.")+id,StateRole::Algebraic);
    auto formulas=mathematics(immersed);
    std::vector<std::string> members;
    for (const auto& formula:formulas) members.push_back(formula.id);
    const auto method=std::string("Immersed.")+FDM::toString(immersed.algorithm);
    ExecutionScope group;
    group.kind=ExecutionKind::Sequence;group.id="immersed.correction";group.order=70;
    for (auto formula:formulas) {
        const std::string target=formula.id=="ibm.momentum"?"rhoU":formula.id=="ibm.energy"?"rhoE"
            :formula.id=="ibm.translation"?"U_s":formula.id=="ibm.rotation"?"omega_s"
            :formula.id=="ibm.incompressibility"?"p":lambda;
        const auto kind=formula.id=="ibm.incompressibility"?TargetKind::Working
            :surface && target==lambda?TargetKind::Workspace:TargetKind::Physical;
        ExecutionScope call;call.kind=ExecutionKind::EquationCall;
        call.step={formula.id,{target,kind}};
        group.children.push_back(std::move(call));
        system.bindNumerics({formula.id,method,members});
        system.addEquation(std::move(formula));
    }
    system.addExecution(std::move(group));
    if (multiplier) system.addConstraint({"immersed.velocity","immersed velocity relation",
        immersed.enforcement==FDM::IBMEnforcement::ExplicitIBM
            ?"lagged multiplier response for the next physical step; not exact current-step no-slip":"J U = U_b / rigid affine constraint",lambda});
}
} // namespace SF::IBM::SystemContribution
