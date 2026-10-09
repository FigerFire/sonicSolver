#include "models/ibm/SF_ibmSystemContribution.h"
#include "models/ibm/descriptor/SF_algorithmDescriptor.h"
#include "models/turbulence/SF_turbulenceSystemContribution.h"
#include "solver/system/SF_immersedMethods.h"
#include "solver/system/SF_systemBuilder.h"
#include "solver/run/SF_planExecutor.h"
#include "core/system/SF_operationIds.h"
#include <algorithm>
#include <memory>
#include <stdexcept>
using namespace SF::System;
namespace {
void require(bool value,const char* message) {if (!value) throw std::runtime_error(message);}
SF::FDM::IBMForcingConfig configuration(SF::FDM::IBMForcingAlgorithm algorithm) {
    using namespace SF::FDM;
    IBMForcingConfig c;c.algorithm=algorithm;
    switch (algorithm) {
    case IBMForcingAlgorithm::PeskinOriginal:
        c.enforcement=IBMEnforcement::ExplicitIBM;c.constraintSupport=IBMConstraintSupport::Surface;break;
    case IBMForcingAlgorithm::DFMExplicitSelfPropelled:
        c.enforcement=IBMEnforcement::ExplicitIBM;c.solidModel=IBMSolidModel::SelfPropelledRigid;break;
    case IBMForcingAlgorithm::DFMFractionalStepSelfPropelled:
        c.enforcement=IBMEnforcement::FractionalDLM;c.solidModel=IBMSolidModel::SelfPropelledRigid;break;
    case IBMForcingAlgorithm::DFMFractionalStepPrescribed: c.enforcement=IBMEnforcement::FractionalDLM;break;
    case IBMForcingAlgorithm::VelocityForcingFTS:
        c.enforcement=IBMEnforcement::VelocityForcing;c.constraintSupport=IBMConstraintSupport::Surface;break;
    case IBMForcingAlgorithm::VelocityForcingBP: c.enforcement=IBMEnforcement::BrinkmanPenalty;break;
    case IBMForcingAlgorithm::DFMImplicitPrescribed:case IBMForcingAlgorithm::DFMAugmentedLagrangian:
        c.enforcement=IBMEnforcement::MonolithicKKT;c.constraintSupport=IBMConstraintSupport::Surface;break;
    case IBMForcingAlgorithm::DFMImplicitSelfPropelled:
        c.enforcement=IBMEnforcement::MonolithicKKT;c.constraintSupport=IBMConstraintSupport::Surface;
        c.solidModel=IBMSolidModel::SelfPropelledRigid;break;
    }
    c.representation=c.constraintSupport==IBMConstraintSupport::Surface?IBMRepresentation::DiffuseKernel:IBMRepresentation::EulerianMask;
    if (c.solidModel==IBMSolidModel::SelfPropelledRigid) c.motion=IBMSolidMotion::CoupledRigid;
    if (c.enforcement==IBMEnforcement::MonolithicKKT) c.coupling=IBMConstraintCoupling::MonolithicKKT;
    c.augmentationCoefficient=100;c.penaltyCoefficient=1000;
    c.constraintSolverMaxIterations=1000;c.constraintSolverRelativeTolerance=1e-8;
    return c;
}
struct BoundSystem : SF::FDM::IImmersedSystem {
    SF::FDM::ImmersedAlgorithmDescriptor descriptor;
    SF::FDM::ImmersedMethodSelection selection;
    SF::FDM::ImmersedMethodCapabilities capability;
    std::vector<SF::FDM::ImmersedFluidPort> ports;
    std::vector<SF::Vector3> force,lagged;
    SF::Vector3 linear,angular;
    explicit BoundSystem(SF::FDM::ImmersedAlgorithmDescriptor d):descriptor(std::move(d)) {
        selection.enforcement=descriptor.enforcement;selection.support=descriptor.support;selection.solid=descriptor.solid;selection.representation=descriptor.representation;
    }
    const SF::FDM::ImmersedMethodSelection& methodSelection() const override {return selection;}
    const SF::FDM::ImmersedMethodCapabilities& capabilities() const override {return capability;}
    const SF::FDM::ImmersedAlgorithmDescriptor& algorithmDescriptor() const override {return descriptor;}
    const std::vector<SF::FDM::ImmersedFluidPort>& fluidPorts() const override {return ports;}
    SF::FDM::ImmersedStorageView storageView(SF::FDM::ImmersedStorageKind kind) const override {
        using K=SF::FDM::ImmersedStorageKind;
        switch(kind) {
        case K::ForceDensity:return {&force,nullptr};
        case K::LaggedForceDensity:return {&lagged,nullptr};
        case K::SolidTranslation:return {nullptr,&linear};
        case K::SolidRotation:return {nullptr,&angular};
        }
        return {};
    }
};
}
int main() {
    using A=SF::FDM::IBMForcingAlgorithm;
    SF::FDM::SolverConfig config;config.numerics.timeRecipe=SF::FDM::builtInTimeRecipe(SF::FDM::TimeRecipeId::ClassicalRK4);
    config.numerics.recipes.time=config.numerics.timeRecipe;
    BuildRequest request;request.singleFluidPreset=SingleFluidPresetSpec{false};
    request.composition.stateDeclared=true;request.composition.solutionVariables={"rho","rhoU","rhoE"};
    std::vector<SF::FDM::IBMForcingConfig> variants;
    for (auto algorithm:{A::PeskinOriginal,A::DFMExplicitSelfPropelled,A::DFMFractionalStepSelfPropelled,
            A::DFMFractionalStepPrescribed,A::VelocityForcingFTS,A::VelocityForcingBP,
            A::DFMImplicitPrescribed,A::DFMImplicitSelfPropelled,A::DFMAugmentedLagrangian})
        variants.push_back(configuration(algorithm));
    for (auto algorithm:{A::PeskinOriginal,A::VelocityForcingFTS}) {
        auto linear=configuration(algorithm);
        linear.surfaceNormalization=SF::FDM::IBMSurfaceNormalization::LinearReproducing;
        variants.push_back(linear);
    }
    for (auto algorithm:{A::DFMExplicitSelfPropelled,A::DFMFractionalStepSelfPropelled}) {
        auto rotated=configuration(algorithm);rotated.rigidMotionMode=SF::FDM::IBMRigidMotionMode::Rotate;
        variants.push_back(rotated);
    }
    auto coupled=configuration(A::DFMImplicitPrescribed);
    coupled.solidModel=SF::FDM::IBMSolidModel::CoupledRigid;coupled.motion=SF::FDM::IBMSolidMotion::CoupledRigid;
    variants.push_back(coupled);
    for (const auto& configuration:variants) {
        const auto algorithm=configuration.algorithm;
        const auto descriptor=SF::IBM::Descriptor::variational(configuration);
        SystemContribution contribution;SF::IBM::SystemContribution::contribute(contribution,descriptor);
        require(contribution.legacyEquations.empty() && contribution.legacyExecution.empty()
            && contribution.policies.empty() && contribution.transformations.empty(),"IBM retained legacy mathematical/execution authority.");
        const bool kkt=descriptor.monolithic;
        auto selected=request;selected.modelContributions={contribution};
        const auto resolved=build(config,selected);
        require(resolved.solvePlan.sourceProgram.legacyEntries.empty(),"IBM HOW retained legacy entries.");
        require(resolved.runtime.report.status==(kkt?RuntimeStatus::Unsupported:RuntimeStatus::Runnable),
            "IBM capability advertised unimplemented KKT predictor or rejected implemented projection.");
        if (!kkt) {
            auto parallel=selected;parallel.parallel=true;
            const auto distributed=build(config,parallel);
            require(distributed.runtime.report.status==RuntimeStatus::Runnable,
                "Implemented distributed surface/body correction was rejected.");
        }
        const auto operation=kkt?OpIds::IbmKktSolve:OpIds::IbmConstraintProject;
        std::vector<std::string> order;SF::Run::OpRegistry callbacks;
        const auto bind=[&](const auto& self,const SolvePlanNode& node)->void {
            require(!node.legacyAdapter,"Native IBM plan contains a compatibility leaf.");
            if (node.operation==OpIds::IbmKktSolve)
                require(node.kind==PlanNodeKind::BlockSolve,"Coupled IBM plan lost block-solve semantics.");
            if (!node.operation.empty()) callbacks.bind(node.operation,node.provider,[&,op=node.operation]{order.push_back(op);});
            for (const auto& child:node.children) self(self,child);
        };
        bind(bind,resolved.solvePlan.root);SF::Run::PlanExecutor::execute(resolved.solvePlan,callbacks);
        require(std::count(order.begin(),order.end(),operation)==1
            && std::count(order.begin(),order.end(),OpIds::ExplicitStageExecute)==4
            && order[7]==operation && order[8]==OpIds::FlowStepCommit,"IBM moved inside RK or after commit.");
        BoundSystem bound(descriptor);validateImmersedBindings(resolved.solvePlan,&bound,false);
        require(bound.storageView(SF::FDM::ImmersedStorageKind::ForceDensity).array==&bound.force,
            "IBM alias copied physical state.");
        auto altered=resolved.solvePlan;
        auto native=std::find_if(altered.compiledProgram.steps.begin(),altered.compiledProgram.steps.end(),
            [](const auto& call) {return call.backendProvider=="ibm.constraint";});
        native->backendOperation=OpIds::FlowStepCommit;
        bool rejected=false;
        try{validateImmersedBindings(altered,&bound,false);}catch(const std::runtime_error&){rejected=true;}
        require(rejected,"Changed numerical operation bypassed frozen IBM binding.");
        if (descriptor.solid==SF::FDM::IBMSolidModel::SelfPropelledRigid && !kkt) {
            bound.descriptor.rigidMotionMode=descriptor.rigidMotionMode==SF::FDM::IBMRigidMotionMode::Rotate
                ?SF::FDM::IBMRigidMotionMode::Motivation:SF::FDM::IBMRigidMotionMode::Rotate;
            rejected=false;
            try{validateImmersedBindings(resolved.solvePlan,&bound,false);}catch(const std::runtime_error&){rejected=true;}
            require(rejected,"Changed rigid DOFs bypassed frozen generalized-mass relation.");
            bound.descriptor.rigidMotionMode=descriptor.rigidMotionMode;
        }
        if (descriptor.support==SF::FDM::IBMConstraintSupport::Surface) {
            const auto original=bound.descriptor.surfaceNormalization;
            bound.descriptor.surfaceNormalization=original==SF::FDM::IBMSurfaceNormalization::PartitionOfUnity
                ?SF::FDM::IBMSurfaceNormalization::LinearReproducing:SF::FDM::IBMSurfaceNormalization::PartitionOfUnity;
            rejected=false;try{validateImmersedBindings(resolved.solvePlan,&bound,false);}catch(const std::runtime_error&){rejected=true;}
            require(rejected,"Runtime surface transfer bypassed frozen WHICH selection.");
            bound.descriptor.surfaceNormalization=original;
        }
        bound.descriptor.algorithm=algorithm==A::PeskinOriginal?A::VelocityForcingFTS:A::PeskinOriginal;
        rejected=false;try{validateImmersedBindings(resolved.solvePlan,&bound,false);}catch(const std::runtime_error&){rejected=true;}
        require(rejected,"Changed runtime IBM implementation bypassed frozen binding.");
        auto user=request;user.userContributions={contribution};
        const auto equivalent=build(config,user);
        require(equivalent.solvePlan.compiledProgram.steps.size()==resolved.solvePlan.compiledProgram.steps.size(),"User/builtin IBM compilation split.");
        for (std::size_t i=0;i<resolved.solvePlan.compiledProgram.steps.size();++i)
            require(equivalent.solvePlan.compiledProgram.steps[i].backendProvider==resolved.solvePlan.compiledProgram.steps[i].backendProvider
                && equivalent.solvePlan.compiledProgram.steps[i].backendOperation==resolved.solvePlan.compiledProgram.steps[i].backendOperation,"Provenance changed IBM execution.");
        const auto fails=[&](SystemContribution bad) {
            auto r=request;r.modelContributions={std::move(bad)};bool failed=false;
            try{(void)build(config,r);}catch(const std::runtime_error&){failed=true;}
            require(failed,"Invalid IBM mathematical/STATE/fusion contract was accepted.");
        };
        auto missing=contribution;missing.execution[0].children.pop_back();fails(missing);
        auto reordered=contribution;std::swap(reordered.execution[0].children[0],reordered.execution[0].children[1]);fails(reordered);
        auto changed=contribution;changed.registeredEquations[0].rhs=FormulaExpr::constantValue(7);fails(changed);
        auto badDescriptor=contribution;badDescriptor.immersed->enforcement=SF::FDM::IBMEnforcement::GhostCell;fails(badDescriptor);
        auto badStorage=contribution;badStorage.states[0].storageKey="newForceCopy";fails(badStorage);
        auto wrongTarget=contribution;wrongTarget.execution[0].children[0].step.target.symbol="rhoE";fails(wrongTarget);
        if (!kkt) {
            auto beforePredictor=request;
            beforePredictor.modelContributions={contribution};
            beforePredictor.authoredExecution=resolved.solvePlan.sourceProgram;
            beforePredictor.authoredNumerics=resolved.numericalSelection.bindings;
            auto& declared=beforePredictor.authoredExecution->root.children;
            const auto correction=std::find_if(declared.begin(),declared.end(),[](const auto& node){return node.id=="immersed.correction";});
            std::iter_swap(declared.begin(),correction);
            bool rejectedPlacement=false;
            try {(void)build(config,beforePredictor);} catch (const std::runtime_error& error) {
                rejectedPlacement=std::string(error.what()).find("placement violation")!=std::string::npos;
            }
            require(rejectedPlacement,"Pre-predictor IBM HOW escaped its declared placement.");
        }
    }
    SystemContribution ghost;SF::IBM::SystemContribution::contribute(ghost,SF::IBM::Descriptor::ghostCell());
    require(ghost.registeredEquations.empty() && ghost.execution.empty() && ghost.boundaryClosures.size()==1,"Ghost fabricated an equation/forcing step.");
    auto withGhost=request;withGhost.modelContributions={ghost};const auto compiledGhost=build(config,withGhost);
    require(compiledGhost.executableSystem.boundaryClosures.size()==1 && compiledGhost.runtime.report.status==RuntimeStatus::Runnable,"Ghost boundary contract not compiled.");
    auto badGhost=withGhost;badGhost.modelContributions[0].boundaryClosures[0].order.pop_back();bool failed=false;
    try{(void)build(config,badGhost);}catch(const std::runtime_error&){failed=true;}require(failed,"Changed Ghost stage order accepted.");
    RawEquationSystem eulerian;eulerian.constraints.push_back({"C_SHARED_PRESSURE","shared pressure","",{}});
    eulerian.boundaryClosures=ghost.boundaryClosures;failed=false;
    try{validateImmersedComposition(eulerian);}catch(const std::runtime_error&){failed=true;}require(failed,"Eulerian phase Ghost was falsely runnable.");
    for (const auto* model:{"kEpsilon","kOmegaSST","Smagorinsky"}) {
        SystemContribution turbulence;SF::Turbulence::contribute(turbulence,{model,{},false});
        for (bool boundary:{false,true}) {
            SystemContribution immersed=ghost;
            if (!boundary) {
                immersed={};SF::IBM::SystemContribution::contribute(immersed,
                    SF::IBM::Descriptor::variational(configuration(A::VelocityForcingBP)));
            }
            auto combined=request;combined.modelContributions={turbulence,immersed};
            const auto system=build(config,combined);
            require(system.rawSystem.registry.contains("mu_t") || system.rawSystem.registry.contains("k"), "Composition lost turbulence mathematics.");
            require(!system.solvePlan.compiledProgram.steps.empty(),"Combination has no native HOW.");
            require(system.runtime.report.reason.find("turbulence.boundary.immersed")!=std::string::npos,"Missing concrete immersed boundary capability was not explained.");
            require(system.runtime.report.status==RuntimeStatus::Unsupported,
                "Turbulence/IBM coupling lacks stage and wall contracts but was advertised as runnable.");
        }
    }
}
