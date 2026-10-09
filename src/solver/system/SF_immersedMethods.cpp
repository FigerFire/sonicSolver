/// @brief IBM provider validates supported AST, binds original storage and numerical operation.
#include "SF_immersedMethods.h"
#include "core/system/SF_operationIds.h"
#include <algorithm>
#include <memory>
#include <stdexcept>
#include <utility>
namespace SF::System {
namespace {
std::vector<Equation> supportedMathematics(const FDM::ImmersedAlgorithmDescriptor& selection) {

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
std::vector<Equation> equations;
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
FDM::ImmersedAlgorithmDescriptor supportedSelection(FDM::IBMForcingAlgorithm algorithm) {
    using A=FDM::IBMForcingAlgorithm;
    FDM::ImmersedAlgorithmDescriptor d;d.algorithm=algorithm;
    switch (algorithm) {
    case A::PeskinOriginal:
        d.enforcement=FDM::IBMEnforcement::ExplicitIBM;d.support=FDM::IBMConstraintSupport::Surface;break;
    case A::DFMExplicitSelfPropelled:
        d.enforcement=FDM::IBMEnforcement::ExplicitIBM;d.solid=FDM::IBMSolidModel::SelfPropelledRigid;break;
    case A::DFMFractionalStepSelfPropelled:
        d.enforcement=FDM::IBMEnforcement::FractionalDLM;d.solid=FDM::IBMSolidModel::SelfPropelledRigid;break;
    case A::DFMFractionalStepPrescribed:
        d.enforcement=FDM::IBMEnforcement::FractionalDLM;break;
    case A::VelocityForcingFTS:
        d.enforcement=FDM::IBMEnforcement::VelocityForcing;d.support=FDM::IBMConstraintSupport::Surface;break;
    case A::VelocityForcingBP:
        d.enforcement=FDM::IBMEnforcement::BrinkmanPenalty;break;
    case A::DFMImplicitPrescribed: case A::DFMAugmentedLagrangian:
        d.enforcement=FDM::IBMEnforcement::MonolithicKKT;d.support=FDM::IBMConstraintSupport::Surface;break;
    case A::DFMImplicitSelfPropelled:
        d.enforcement=FDM::IBMEnforcement::MonolithicKKT;d.support=FDM::IBMConstraintSupport::Surface;
        d.solid=FDM::IBMSolidModel::SelfPropelledRigid;break;
    }
    d.representation=d.support==FDM::IBMConstraintSupport::Surface?FDM::IBMRepresentation::DiffuseKernel:FDM::IBMRepresentation::EulerianMask;
    return d;
}
class ImmersedMethod final : public IProvider {
public:
    explicit ImmersedMethod(FDM::IBMForcingAlgorithm algorithm,
                            FDM::IBMSurfaceNormalization normalization=FDM::IBMSurfaceNormalization::PartitionOfUnity)
        : selection_(supportedSelection(algorithm)),id_(std::string("Immersed.")+FDM::toString(algorithm)
            +(normalization==FDM::IBMSurfaceNormalization::LinearReproducing?".linearReproducing":"")) {
        selection_.surfaceNormalization=normalization;
    }
    std::string_view id() const override {return id_;}
    std::string_view runtimeProvider() const override {return "ibm.constraint";}
    CompiledEquationCall compile(const ExecutableEquationSystem& system,const EquationCall& call,
                                const NumericalBinding& binding) const override {
        if (!system.immersed) throw std::runtime_error("Missing frozen immersed descriptor for selected provider.");
        const auto selection=*system.immersed;
        const bool supportedSolid=selection.solid==selection_.solid
            || (selection_.algorithm==FDM::IBMForcingAlgorithm::DFMImplicitPrescribed && selection.solid==FDM::IBMSolidModel::CoupledRigid);
        if (selection.algorithm!=selection_.algorithm || selection.enforcement!=selection_.enforcement
            || selection.support!=selection_.support || selection.representation!=selection_.representation
            || selection.surfaceNormalization!=selection_.surfaceNormalization || !supportedSolid)
            throw std::runtime_error("Selected immersed provider does not support resolved descriptor properties: "+id_);
        auto expected=supportedMathematics(selection);
        const auto matches=[&](const auto& formulas) {
            return std::all_of(formulas.begin(),formulas.end(),[&](const auto& math) {
                return system.registry.contains(math.id)
                    && canonicalFormula(system.registry.at(math.id))==canonicalFormula(math);
            });
        };
        if (!matches(expected)) throw std::runtime_error("Unsupported immersed mathematics for resolved descriptor: "+id_);
        std::vector<std::string> ids;std::vector<CompiledMathRef> members;
        const bool surface=selection.support==FDM::IBMConstraintSupport::Surface;
        const std::string lambda=surface?"Lambda_s":"lambda_b";
        for (const auto& math:expected) {
            ids.push_back(math.id);
            if (canonicalFormula(system.registry.at(math.id))!=canonicalFormula(math))
                throw std::runtime_error("Unsupported immersed mathematics: "+math.id+" for selected "+id_);
            members.push_back({math.id,math.id=="ibm.momentum"?"rhoU":math.id=="ibm.energy"?"rhoE"
                :math.id=="ibm.translation"?"U_s":math.id=="ibm.rotation"?"omega_s"
                :math.id=="ibm.incompressibility"?"p":lambda});
        }
        if (binding.inputs!=ids)
            throw std::runtime_error("Immersed method requires its complete ordered impulse/work/constraint block.");
        const auto found=std::find_if(members.begin(),members.end(),[&](const auto& m) {return m.equation==call.equation;});
        const auto kind=call.equation=="ibm.incompressibility"?TargetKind::Working
            :surface && call.target.symbol==lambda?TargetKind::Workspace:TargetKind::Physical;
        if (found==members.end() || call.target.symbol!=found->target || call.target.kind!=kind)
            throw std::runtime_error("Immersed method has a wrong semantic target: "+call.equation);
        for (const auto& packed:std::vector<std::pair<std::string,int>>{{"rho",0},{"rhoU",1},{"rhoE",4}}) {
            const auto& state=system.state.at(packed.first);
            if (state.storageBinding!=StorageBinding::PackedDistributed || state.storageKey!="conservative"
                || state.componentOffset!=packed.second || state.components!=(packed.first=="rhoU"?3:1))
                throw std::runtime_error("Unsupported immersed fluid storage: "+packed.first+"; current port requires conservative state.");
        }
        for (const auto& name:std::vector<std::string>{"ibm.forceDensity"}) {
            const auto& state=system.state.at(name);
            if (state.storageBinding!=StorageBinding::SpecializedExecutor || state.storageKey!=name || state.components!=3
                || state.role!=StateRole::Derived || state.location!=VariableLocation::EulerianCell
                || state.ownership!=OwnershipKind::EulerianGlobalDof)
                throw std::runtime_error("Invalid immersed original-storage binding: "+name);
        }
        for (const auto& member:members) {
            if (member.target=="rhoU" || member.target=="rhoE" || member.target=="p") continue;
            const auto& state=system.state.at(member.target);
            const auto storage=member.target==lambda ? (surface?"ibm.surfaceMultiplier":"ibm.forceDensity")
                :"ibm."+member.target;
            if (state.components!=3 || state.storageKey!=storage
                || state.role!=(member.target==lambda?StateRole::Multiplier:StateRole::Algebraic)
                || state.location!=(member.target==lambda?(surface?VariableLocation::SurfaceConstraint:VariableLocation::BodyConstraint):VariableLocation::SolidGlobal)
                || state.ownership!=(member.target==lambda?OwnershipKind::ConstraintGlobalDof:OwnershipKind::SolidGlobalDof)
                || state.storageBinding!=(surface && member.target==lambda?StorageBinding::TransientWorkspace:StorageBinding::SpecializedExecutor))
                throw std::runtime_error("Invalid immersed multiplier/solid storage: "+member.target);
        }
        if (selection.enforcement==FDM::IBMEnforcement::ExplicitIBM) {
            const auto& state=system.state.at("ibm.laggedForceDensity");
            if (state.storageKey!="ibm.laggedForceDensity" || state.storageBinding!=StorageBinding::SpecializedExecutor || state.components!=3)
                throw std::runtime_error("Missing original lagged force storage.");
        }
        CompiledEquationCall result;result.source=call;result.target.symbol=call.target.symbol;result.target.kind=call.target.kind;
        result.target.viewOwner=StateViewOwner::NumericalProvider;
        if (kind==TargetKind::Workspace) {
            result.target.workspace="ibm.surfaceMultiplier";
            result.target.workspaceLocation=VariableLocation::SurfaceConstraint;
            result.target.workspaceOwnership=OwnershipKind::ConstraintGlobalDof;
            result.target.resources={{lambda,"ibm.surfaceMultiplier",0,3,ResourceAccessMode::Write,false,SynchronizationRequirement::WriteOwned}};
        } else if (kind==TargetKind::Working) result.target.workspace="ibm.pCorrection";
        result.equationMethod=id_;result.calls={{call.equation,call.target.symbol}};
        result.capabilities={"velocity.immersed"};
        result.stateUses={{"rho"},{"rhoU"},{"rhoE"},{"U"}};
        result.stateEffects={{"rhoU",false,true},{"rhoE",false,true}};
        result.reads={"rho","rhoU","rhoE","U","geometry","targetTime","dt","U_b"};
        if (selection.enforcement==FDM::IBMEnforcement::BrinkmanPenalty) {
            result.reads.push_back("ibm.mask");result.reads.push_back("ibm.penaltyCoefficient");
        }
        if (selection.solid!=FDM::IBMSolidModel::Prescribed) {
            result.reads.insert(result.reads.end(),{"U_s","omega_s","U_deform","solid.externalForce","solid.externalTorque"});
        }
        if (selection.enforcement==FDM::IBMEnforcement::ExplicitIBM) {
            result.reads.push_back("ibm.laggedForceDensity");
            result.stateUses.push_back({"ibm.laggedForceDensity",StateVersion::Lagged});
        }
        // One fused numerical call writes the complete mathematical block.
        result.writes={"ibm.forceDensity"};
        for (const auto& member:members) result.writes.push_back(member.target);
        if (selection.enforcement==FDM::IBMEnforcement::ExplicitIBM) result.writes.push_back("ibm.laggedForceDensity");
        result.fusionKey=id_+"/completeBlock";result.fusionMembers=members;
        const auto operation=selection.enforcement==FDM::IBMEnforcement::MonolithicKKT?OpIds::IbmKktSolve:OpIds::IbmConstraintProject;
        if (selection.enforcement==FDM::IBMEnforcement::MonolithicKKT) {
            result.fragment.kind=PlanNodeKind::BlockSolve;
            result.fragment.id="immersed.block";result.fragment.name="coupled immersed block solve";
            result.fragment.operation=operation;result.fragment.equationCalls=result.calls;
        } else result.backendOperation=operation;
        result.operations.push_back({operation,"bound immersed numerical block",OperationStage::CorrectionCommit,
            {OperationCapability::ImmersedConstraint},{OriginKind::Generated,id_}});
        result.requirements={"complete post-predictor impulse/work block","targetTime = physical time + dt",
            "original multiplier/solid storage; no physical copy","owner COPY for multiplier, SUM for load"};
        result.providerContract=std::make_shared<const CompiledImmersedContract>(CompiledImmersedContract{
            selection.algorithm,selection.enforcement,selection.support,selection.solid,selection.representation,selection.rigidMotionMode,selection.surfaceNormalization,kind,call.equation,call.target.symbol,members});
        return result;
    }
private:
    FDM::ImmersedAlgorithmDescriptor selection_;
    std::string id_;
};
} // namespace
void addImmersedMethods(ProviderRegistry& registry) {
    using A=FDM::IBMForcingAlgorithm;
    static const std::vector<ImmersedMethod> methods={ImmersedMethod(A::PeskinOriginal),ImmersedMethod(A::DFMExplicitSelfPropelled),
        ImmersedMethod(A::DFMFractionalStepSelfPropelled),ImmersedMethod(A::DFMFractionalStepPrescribed),
        ImmersedMethod(A::VelocityForcingFTS),ImmersedMethod(A::VelocityForcingBP),
        ImmersedMethod(A::DFMImplicitPrescribed),ImmersedMethod(A::DFMImplicitSelfPropelled),ImmersedMethod(A::DFMAugmentedLagrangian)};
    for (const auto& method:methods) registry.add(method);
    using N=FDM::IBMSurfaceNormalization;
    static const std::vector<ImmersedMethod> linearMethods={
        ImmersedMethod(A::PeskinOriginal,N::LinearReproducing),ImmersedMethod(A::VelocityForcingFTS,N::LinearReproducing),
        ImmersedMethod(A::DFMImplicitPrescribed,N::LinearReproducing),ImmersedMethod(A::DFMImplicitSelfPropelled,N::LinearReproducing),
        ImmersedMethod(A::DFMAugmentedLagrangian,N::LinearReproducing)};
    for (const auto& method:linearMethods) registry.add(method);
}
void validateImmersedComposition(const RawEquationSystem& raw) {
    const bool ghost=std::any_of(raw.boundaryClosures.begin(),raw.boundaryClosures.end(),
        [](const auto& closure) {return closure.provider=="ibm.boundary";});
    const bool correction=raw.registry.contains("ibm.momentum");
    if (!ghost && !correction) return;
    if (!raw.state.contains("rhoU") && !ghost)
        throw std::runtime_error("Unsupported immersed fluid port: participating STATE has no implemented conservative momentum target rhoU; phase-wise fluid-port assembly is unavailable.");
    for (const auto& closure:raw.boundaryClosures) {
        if (closure.provider!="ibm.boundary") continue;
        if (closure.reads!=std::vector<std::string>{"conservative","geometry","classification"}
            || closure.writes!=std::vector<std::string>{"conservative.ghost"}
            || !closure.atEverySpatialEvaluation
            || closure.order!=std::vector<std::string>{"physical boundary","halo COPY","ghost/ILW reconstruction","ghost publication","halo COPY","spatial operator"})
            throw std::runtime_error("Unsupported: selected Ghost implementation cannot realize changed boundary contract.");
        if (!raw.state.contains("rhoU") || raw.state.at("rhoU").storageBinding!=StorageBinding::PackedDistributed)
            throw std::runtime_error("Unsupported: Ghost boundary requires its implemented conservative ghost-state view.");
    }
}
void validateImmersedBindings(const CompiledSolvePlan& plan,const FDM::IImmersedSystem* system,bool boundaryPresent) {
    for (const auto& call:plan.compiledProgram.steps) {
        const auto* contract=std::any_cast<std::shared_ptr<const CompiledImmersedContract>>(&call.providerContract);
        if (!contract) {
            if (call.backendProvider=="ibm.constraint") throw std::runtime_error("Native IBM occurrence lost its immutable provider contract.");
            continue;
        }
        if (!*contract || !system) throw std::runtime_error("Native immersed occurrence has no matching bound IBM system.");
        const auto& c=**contract;const auto& d=system->algorithmDescriptor();
        for (const auto& port:system->fluidPorts()) if (port.phaseIndex>=0
            || port.momentumName!=c.fluid.momentumName || port.energyName!=c.fluid.energyName
            || port.density!=c.fluid.density || port.velocity!=c.fluid.velocity
            || port.momentumTarget!=c.fluid.momentumTarget || port.energyTarget!=c.fluid.energyTarget
            || !port.coupleMechanicalWork || !port.constrainVelocity)
            throw std::runtime_error("Unsupported bound immersed fluid port: no matching frozen conservative velocity/momentum/energy storage contract.");
        const auto& selection=system->methodSelection();
        const bool sameMembers=call.fusionMembers.size()==c.members.size()
            && std::equal(c.members.begin(),c.members.end(),call.fusionMembers.begin(),
                [](const auto& a,const auto& b) {return a.equation==b.equation && a.target==b.target;});
        const auto force=system->storageView(FDM::ImmersedStorageKind::ForceDensity);
        if (!force.array) throw std::runtime_error("IBM force STATE has no original storage alias.");
        if (c.enforcement==FDM::IBMEnforcement::ExplicitIBM
            && !system->storageView(FDM::ImmersedStorageKind::LaggedForceDensity).array)
            throw std::runtime_error("IBM lagged force STATE has no original storage alias.");
        if ((c.solid==FDM::IBMSolidModel::SelfPropelledRigid || c.solid==FDM::IBMSolidModel::CoupledRigid)
            && (!system->storageView(FDM::ImmersedStorageKind::SolidTranslation).value
                || !system->storageView(FDM::ImmersedStorageKind::SolidRotation).value))
            throw std::runtime_error("IBM solid STATE has no original velocity alias.");
        if (d.surfaceNormalization!=c.surfaceNormalization || d.algorithm!=c.algorithm
            || (c.solid==FDM::IBMSolidModel::SelfPropelledRigid && c.enforcement!=FDM::IBMEnforcement::MonolithicKKT
                && d.rigidMotionMode!=c.rigidMotionMode)
            || selection.enforcement!=c.enforcement || selection.support!=c.support
            || selection.solid!=c.solid || selection.representation!=c.representation
            || call.source.target.kind!=c.kind || call.target.kind!=c.kind || call.source.target.symbol!=c.target
            || call.equationMethod!=std::string("Immersed.")+FDM::toString(c.algorithm)
                +(c.surfaceNormalization==FDM::IBMSurfaceNormalization::LinearReproducing?".linearReproducing":"")
            || call.backendProvider!="ibm.constraint" || call.source.equation!=c.equation || call.target.symbol!=c.target
            || !sameMembers || (c.enforcement==FDM::IBMEnforcement::MonolithicKKT
                ?(!call.backendOperation.empty() || call.fragment.kind!=PlanNodeKind::BlockSolve
                    || call.fragment.operation!=OpIds::IbmKktSolve || !call.fragment.children.empty())
                :call.backendOperation!=OpIds::IbmConstraintProject))
            throw std::runtime_error("Bound IBM implementation differs from frozen WHAT/STATE/HOW/WHICH contract.");
    }
    if (boundaryPresent && system && system->methodSelection().enforcement!=FDM::IBMEnforcement::GhostCell)
        throw std::runtime_error("Compiled immersed ghost boundary has no matching bound Ghost system.");
}
} // namespace SF::System
