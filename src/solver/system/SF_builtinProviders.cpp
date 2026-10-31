#include "core/interfaces/SF_transportModel.h"
#include "SF_fusedSpatial.h"
#include "models/turbulence/SF_transportMathematics.h"
#include "SF_thermalMethod.h"
/// Selected numerical implementations compile mathematics; HOW owns the process.
#include "SF_methodObjects.h"
#include "SF_immersedMethods.h"
#include "SF_levelSetMethods.h"
#include "SF_multiphaseMethods.h"
#include "SF_eulerianCoupling.h"
#include "SF_pressureCoupling.h"
#include "SF_scalarMethod.h"
#include "SF_formulaCompiler.h"
#include "core/system/SF_operationIds.h"
#include <algorithm>
#include <stdexcept>
#include <utility>
namespace SF::System {
namespace {

void bindSourceInputs(const ExecutableEquationSystem& system,const EquationCall& source,
        const Equation& equation,CompiledEquationCall& call) {
    const auto visit=[&](const auto& self,const FormulaExpr& e)->void {
        if(e.kind()==FormulaExpr::Kind::Operator && e.operation()==Math::Op::Source && e.arguments().size()==1) {
            for(const auto& id:inspectFormula(e.arguments()[0]).symbols) {
                if(!system.state.contains(id) || system.state.at(id).constantValue)continue;
                const auto selected=std::find_if(source.inputs.begin(),source.inputs.end(),[&](const auto& read){return read.symbol==id;});
                if(selected==source.inputs.end())throw std::runtime_error("Field source requires explicit STATE input: step="+source.occurrence+" equation="+equation.id+" input="+id);
                if(selected->kind!=TargetKind::Physical || selected->halo || !selected->snapshot.empty() || selected->perIteration)
                    throw std::runtime_error("Unsupported field source input qualifier: "+id);
                const auto previous=std::find_if(call.stateUses.begin(),call.stateUses.end(),[&](const auto& read){return read.symbol==id;});
                if(previous==call.stateUses.end())call.stateUses.push_back(*selected);
                else if(previous->version!=selected->version)throw std::runtime_error("Conflicting field source input versions: "+id);
            }
            return;
        }
        for(const auto& child:e.arguments())self(self,child);
    };
    visit(visit,equation.rhs);
}

SolvePlanNode methodLeaf(PlanNodeKind kind,std::string id,
                         std::string name,std::string operation,
                         std::string_view provider = {}) {
    SolvePlanNode leaf;
    leaf.kind=kind;
    leaf.id=std::move(id);
    leaf.name=std::move(name);
    leaf.operation=std::move(operation);
    leaf.provider=provider;
    return leaf;
}

SolvePlanNode explicitFragment(const CompiledEquationCall& body,
                               const CompiledTimeRecipe& recipe,
                               const std::vector<SolvePlanNode>& prefix) {
    if (!body.temporalResidual || body.calls.empty()
        || body.backendOperation.empty())
        throw std::runtime_error("TemporalMethod requires a compiled residual body and backend.");
    SolvePlanNode root;
    root.kind=PlanNodeKind::Sequence;
    root.id="Explicit.step";
    root.name="explicit time step";
    if (!body.temporalParticipants.empty()) {
        const auto leaf=[&](const std::string& op,const std::string& owner,
                const TemporalParticipant* participant=nullptr) {
            auto node=methodLeaf(PlanNodeKind::Update,op,op,op,owner);
            if (participant) {node.occurrence=participant->identity;node.equationCalls=participant->calls;}
            return node;
        };
        for (const auto& p:body.temporalParticipants) for (const auto& op:p.preparation)
            root.children.push_back(leaf(op,p.provider,&p));
        root.children.push_back(leaf(TemporalOps::Dt,TemporalOps::Provider));
        root.children.insert(root.children.end(),prefix.begin(),prefix.end());
        for (const auto& p:body.temporalParticipants) root.children.push_back(leaf(p.snapshot,p.provider,&p));
        auto stages=methodLeaf(PlanNodeKind::StageLoop,"Explicit.stages","synchronous explicit stages",{});
        stages.repetitions=recipe.stageCount();
        stages.children.push_back(leaf(TemporalOps::Open,TemporalOps::Provider));
        for (const auto& p:body.temporalParticipants) stages.children.push_back(leaf(p.prepareStage,p.provider,&p));
        stages.children.push_back(leaf(TemporalOps::Ready,TemporalOps::Provider));
        for(const auto& p:body.temporalParticipants)if(!p.coefficientPrepare.empty())stages.children.push_back(leaf(p.coefficientPrepare,p.provider,&p));
        for (const auto& p:body.temporalParticipants) stages.children.push_back(leaf(p.rhs,p.provider,&p));
        stages.children.push_back(leaf(TemporalOps::RhsReady,TemporalOps::Provider));
        for (const auto& p:body.temporalParticipants) stages.children.push_back(leaf(p.advance,p.provider,&p));
        stages.children.push_back(leaf(TemporalOps::Close,TemporalOps::Provider));
        root.children.push_back(std::move(stages));return root;
    }
    for (const auto& operation:body.temporalPreparation)
        root.children.push_back(methodLeaf(PlanNodeKind::Update,
            "Explicit.prepare","prepare physical step",operation,body.backendProvider));
    if (body.temporalStepSize.empty() || body.temporalSnapshot.empty())
        throw std::runtime_error("Temporal provider must declare step-size and snapshot operations.");
    root.children.push_back(methodLeaf(PlanNodeKind::Update,
        "Explicit.dt","compute time step",body.temporalStepSize,body.backendProvider));
    root.children.insert(root.children.end(),prefix.begin(),prefix.end());
    root.children.push_back(methodLeaf(PlanNodeKind::Update,
        "Explicit.begin","begin explicit step",body.temporalSnapshot,body.backendProvider));
    auto stages=methodLeaf(PlanNodeKind::StageLoop,"Explicit.stages",
        std::string(FDM::toString(recipe.id()))+" fused explicit stages",{});
    stages.repetitions=recipe.stageCount();
    auto stage=methodLeaf(PlanNodeKind::Update,"Explicit.stage.execute",
        "execute explicit stage",body.backendOperation,body.backendProvider);
    stage.equationCalls=body.calls;
    stages.children.push_back(std::move(stage));
    root.children.push_back(std::move(stages));
    return root;
}

class ForwardEulerMethod final : public ITemporalMethod {
public:
    FDM::TimeRecipeId id() const override { return FDM::TimeRecipeId::ForwardEuler; }
    CompiledTimeRecipe compile(FDM::TimeRecipe selected) const override {
        return CompiledTimeRecipe::fromMethod(selected,{{{0.0,0.0,1.0}}},1,
            ExplicitStageBackend::ForwardEuler);
    }
    SolvePlanNode compileFragment(const CompiledTimeRecipe& compiled,
                                  const CompiledEquationCall& body,
                                  const std::vector<SolvePlanNode>& prefix) const override {
        return explicitFragment(body,compiled,prefix);
    }
};

class SSPRK3Method final : public ITemporalMethod {
public:
    FDM::TimeRecipeId id() const override { return FDM::TimeRecipeId::SSPRK3; }
    CompiledTimeRecipe compile(FDM::TimeRecipe selected) const override {
        return CompiledTimeRecipe::fromMethod(selected,
            {{{0.0,0.0,1.0},{1.0,0.75,0.25},{0.5,1.0/3.0,2.0/3.0}}},3,
            ExplicitStageBackend::SSPRK3);
    }
    SolvePlanNode compileFragment(const CompiledTimeRecipe& compiled,
                                  const CompiledEquationCall& body,
                                  const std::vector<SolvePlanNode>& prefix) const override {
        return explicitFragment(body,compiled,prefix);
    }
};

class ClassicalRK4Method final : public ITemporalMethod {
public:
    FDM::TimeRecipeId id() const override { return FDM::TimeRecipeId::ClassicalRK4; }
    CompiledTimeRecipe compile(FDM::TimeRecipe selected) const override {
        return CompiledTimeRecipe::fromMethod(selected,
            {{{0.0,0.0,0.5},{0.5,0.0,0.5},{0.5,0.0,1.0},{1.0,0.0,0.0}}},
            4,ExplicitStageBackend::ClassicalRK4,
            {{1.0,2.0,2.0,1.0}},6.0);
    }
    SolvePlanNode compileFragment(const CompiledTimeRecipe& compiled,
                                  const CompiledEquationCall& body,
                                  const std::vector<SolvePlanNode>& prefix) const override {
        return explicitFragment(body,compiled,prefix);
    }
};

bool containsDdt(const FormulaExpr& expression, std::string_view target) {
    if (expression.kind()==FormulaExpr::Kind::Operator
        && expression.operation()==Math::Op::Ddt && expression.arguments().size()==1
        && expression.arguments().front().kind()==FormulaExpr::Kind::Symbol
        && expression.arguments().front().name()==target) return true;
    return std::any_of(expression.arguments().begin(),expression.arguments().end(),
        [&](const FormulaExpr& child) { return containsDdt(child,target); });
}

void collectSymbols(const FormulaExpr& expression,std::vector<std::string>& result) {
    for(const auto& symbol:inspectFormula(expression).symbols)
        if(std::find(result.begin(),result.end(),symbol)==result.end())result.push_back(symbol);
}

void requireOutput(const ExecutableEquationSystem&,const EquationCall& step) {
    if (step.equation.empty() || step.target.symbol.empty())
        throw std::runtime_error("EquationCall requires an equation and a STATE target.");
}

CompiledTarget bindTarget(const ExecutableEquationSystem&,const EquationCall& call,
                          std::string workspace = {}) {
    CompiledTarget result;
    result.symbol=call.target.symbol;result.kind=call.target.kind;
    if (call.target.kind!=TargetKind::Physical)
        result.workspace=workspace.empty()
            ? "state."+call.target.symbol+"."+(call.target.kind==TargetKind::Working
                ? "Working" : call.target.kind==TargetKind::Correction ? "Correction" : "Workspace")
            : std::move(workspace);
    return result;
}

/// The existing nonlocal kernel consumes this complete mathematical group once.
class SstClosureMethod final : public IProvider {
public:
    std::string_view id() const override{return "SSTNonlocalClosure";}
    std::string_view runtimeProvider() const override{return ScalarOps::ClosureProvider;}
    CompiledEquationCall compile(const ExecutableEquationSystem& system,const EquationCall& call,const NumericalBinding& binding) const override {
        const auto equations=Turbulence::stageClosureMathematics();
        std::vector<std::string> outputs;for(const auto& e:equations)outputs.push_back(e.id);
        if(!call.stageScoped || call.target.symbol!=outputs.front() || binding.inputs!=outputs || !binding.parameters.empty() || !binding.operators.empty())throw std::runtime_error("SST nonlocal implementation requires its complete declared output group in Stage scope");
        CompiledEquationCall result;result.source=call;result.target.symbol=call.target.symbol;result.stageLocal=true;
        result.equationMethod=id();result.backendOperation=instanceOperation("closure.evaluate",call.occurrence);
        result.temporalPublish=instanceOperation("closure.publish",call.occurrence);
        result.coefficientPrepare=result.backendOperation;result.coefficientPublish=result.temporalPublish;
        result.coefficientProvides={instanceOperation("closure.outputs",call.occurrence)};result.coefficientStates=outputs;
        result.stateUses={{"rho",StateVersion::Stage},{"rhoU",StateVersion::Stage},{"k",StateVersion::Stage},{"omega",StateVersion::Stage}};
        if(call.inputs.size()!=result.stateUses.size() || !std::equal(call.inputs.begin(),call.inputs.end(),result.stateUses.begin(),[](const auto& a,const auto& b){return a.symbol==b.symbol && a.version==b.version && a.kind==TargetKind::Physical && !a.halo && a.snapshot.empty() && !a.perIteration;}))throw std::runtime_error("SST closure requires explicit Stage rho/rhoU/k/omega inputs");
        result.coefficientInputs=result.stateUses;
        for(const auto& e:equations) {
            result.operatorBindings.push_back(e.rhs.name());
            if(canonicalFormula(system.registry.at(e.id))!=canonicalFormula(e))throw std::runtime_error("Unsupported nonlocal SST mathematics: "+e.id);
            system.state.at(e.id);result.calls.push_back({e.id,e.id});result.writes.push_back(e.id);result.stateEffects.push_back({e.id,true,false,true});
        }
        result.providerContract=CompiledStageClosure{{"k","omega"}};result.requiredHaloWidth=1;
        result.physicalResources={{call.occurrence,"model.turbulence",typeid(FDM::ITransportModel),true}};
        for(const auto& op:{result.backendOperation,result.temporalPublish})result.operations.push_back({op,"evaluate complete nonlocal SST outputs",OperationStage::Prepare,{OperationCapability::ScalarExplicit},{OriginKind::Model,"SST nonlocal"}});
        result.requirements={"serial single patch; stationary non-IBM stage port","original SST gradients, wall distance, molecular viscosity and regularization; complete output group evaluated once","not a temporal unknown; authored Stage evaluation and final Current refresh"};return result;
    }
};
CompiledEquationCall compileConservativeResidual(const ExecutableEquationSystem& system,const EquationCall& step,const NumericalBinding& binding) {
        requireOutput(system,step);
        if(!binding.inputs.empty() || !binding.parameters.empty())
            throw std::runtime_error("ConservativeResidual consumes its actual equation and STATE inputs, not auxiliary mathematics or hidden parameters: "+step.occurrence);
        if (step.target.kind!=TargetKind::Physical)
            throw std::runtime_error("Unsupported: conservative residual provider requires a Physical target.");
        if (step.target.symbol.empty())
            throw std::runtime_error("ConservativeResidual requires one equation and one target.");
        CompiledEquationCall result;
        result.source=step;
        result.target=bindTarget(system,step);
        result.equationMethod="ConservativeResidual";
        result.physicalPublication={OpIds::FlowStepCommit,"flow.conservative"};
        result.temporalResidual=true;
        result.temporalPreparation={OpIds::FlowStepPrepare};
        result.temporalStepSize=OpIds::FlowDtCompute;
        result.temporalSnapshot=OpIds::FlowStepBegin;
        result.temporalOwner="flow.conservative";
        if(step.target.symbol=="rho")result.stageProvides={"canonicalMassFlux(Stage)"};
        result.temporalStagePrepare="flow.stage.prepare";result.temporalRhs="flow.stage.rhs";
        result.temporalAdvance="flow.stage.advance";result.temporalPublish=OpIds::FlowStepCommit;
        result.oldTimeWorkspace="explicit.q0";
        result.publishesStageToPhysicalTarget=true;
        if (system.state.contains("mu_t")) {
            result.coefficientUses.push_back({"mu_t",{StateVersion::Frozen,StateVersion::Stage}});
        }
        result.spatialTerms=true;
        result.residualWorkspace="residual";
        for (const char* operation:{OpIds::FlowStepPrepare,OpIds::FlowDtCompute,OpIds::FlowStepBegin,
                OpIds::ExplicitStageExecute,OpIds::FlowStepCommit,
                "flow.stage.prepare","flow.stage.rhs","flow.stage.advance"})
            result.operations.push_back({operation,operation,OperationStage::Prepare,
                {OperationCapability::ConservativeExplicit},{OriginKind::Generated,"ConservativeResidual"}});
        result.backendOperation=OpIds::ExplicitStageExecute;
        result.writes={step.target.symbol};
        result.requirements={"boundary-ready state","conservative residual"};
        const auto& formula=system.registry.at(step.equation);
        const auto& symbol=step.target.symbol;
        if (!containsDdt(formula.lhs,symbol))
            throw std::runtime_error("ConservativeResidual equation has no ddt(target).");
        result.calls.push_back({formula.id,symbol});
        collectSymbols(formula.lhs,result.reads);
        collectSymbols(formula.rhs,result.reads);
        const auto stageSources=[&](const auto& self,const FormulaExpr& expression)->void {
            if (expression.kind()==FormulaExpr::Kind::Operator && expression.operation()==Math::Op::Source && expression.arguments().size()==4
                && expression.arguments()[0].kind()==FormulaExpr::Kind::Symbol
                && (expression.arguments()[0].name()=="linearScalarForce" || expression.arguments()[0].name()=="linearScalarWork")) {
                const auto& symbol=system.state.at(expression.arguments()[1].name());
                if (symbol.components!=1 || symbol.storageBinding!=StorageBinding::NamedDistributed || symbol.readOnly())
                    throw std::runtime_error("Scalar feedback requires scalar NamedDistributed Stage input.");
                result.stateUses.push_back({symbol.id,StateVersion::Stage});
                if (expression.arguments()[0].name()=="linearScalarWork") {
                    result.stateUses.push_back({"rho",StateVersion::Stage});
                    result.stateUses.push_back({"rhoU",StateVersion::Stage});
                }
            }
            for(const auto& child:expression.arguments())self(self,child);
        };
        stageSources(stageSources,formula.rhs);
        bindSourceInputs(system,step,formula,result);
        return result;
    }// Capability matcher for the frozen backend, not an alternate WHAT registry.
Equation supportedRasMathematics(const std::string& second,
                                            const std::string& target) {
    using E=FormulaExpr;
    const auto s=[](const std::string& name) { return Turbulence::mathematicalSymbol(name); };
    const auto mul=[](E a,E b) { return E::multiply(std::move(a),std::move(b)); };
    const auto div=[](E a,E b) { return E::divide(std::move(a),std::move(b)); };
    if ((second!="epsilon" && second!="omega") || (target!="k" && target!=second))
        throw std::runtime_error("Unsupported single-fluid RAS mathematical contract.");
    const bool epsilon=second=="epsilon";
    const std::string model=epsilon ? "kEpsilon" : "kOmegaSST";
    const auto coefficient=[&](const std::string& name) { return s(model+"."+name); };
    const auto production=coefficient("limitedProduction");
    E source;
    E diffusivity;
    if (target=="k") {
        const auto sink=epsilon ? s("epsilon")
            : mul(coefficient("betaStar"),mul(s("k"),s("omega")));
        source=E::subtract(div(production,s("rho")),sink);
    } else if (epsilon) {
        source=E::subtract(div(mul(coefficient("C1"),mul(s("epsilon"),production)),
                                  coefficient("maxRhoK")),
            div(mul(coefficient("C2"),mul(s("epsilon"),s("epsilon"))),coefficient("maxK")));
    } else {
        source=E::subtract(div(mul(coefficient("gamma"),production),coefficient("maxMuT")),
            mul(coefficient("beta"),mul(s("omega"),s("omega"))));
    }
    diffusivity=epsilon ? div(s("mu_t"),coefficient(target=="k" ? "sigmaK" : "sigmaEpsilon"))
        : E::add(coefficient("laminarMu"),mul(coefficient(target=="k" ? "sigmaK" : "sigmaOmega"),s("mu_t")));
    auto diffusion=E::op(Math::Op::Diffusion,{diffusivity,s(target)});
    if (!epsilon && target=="omega")
        diffusion=E::add(std::move(diffusion),coefficient("crossDiffusion"));
    return {target,E::op(Math::Op::Ddt,{s(target)}),
        E::add(std::move(source),div(std::move(diffusion),s("rho"))),
        {OriginKind::Model,model}};
}

class TurbulenceTransportMethod final : public IProvider {
public:
    std::string_view id() const override { return "TurbulenceTransport"; }
    std::string_view runtimeProvider() const override { return "flow.turbulence"; }
    CompiledEquationCall compile(const ExecutableEquationSystem& system,
                                const EquationCall& step,
                                const NumericalBinding& binding) const override {
        if (binding.inputs.size()!=2 || binding.inputs[0]!="k"
            || (binding.inputs[1]!="omega" && binding.inputs[1]!="epsilon"))
            throw std::runtime_error("TurbulenceTransport requires an explicit complete RAS pair.");
        const auto& second=binding.inputs[1];
        const auto model=second=="omega" ? "kOmegaSST" : "kEpsilon";
        if (system.state.contains(second=="omega" ? "epsilon" : "omega"))
            throw std::runtime_error("TurbulenceTransport has incompatible model/state pair.");
        if (step.target.kind!=TargetKind::Physical || step.target.symbol!=step.equation
            || (step.equation!="k" && step.equation!=second))
            throw std::runtime_error("TurbulenceTransport requires each equation's own Physical target.");
        for (const auto& name:binding.inputs) {
            const auto& state=system.state.at(name);
            if (state.role!=StateRole::Transported || state.storageBinding!=StorageBinding::ProviderDistributed
                || state.storageKey!=name || state.components!=1)
                throw std::runtime_error("TurbulenceTransport requires provider-owned transported STATE: "+name);
            if (canonicalFormula(system.registry.at(name))
                !=canonicalFormula(supportedRasMathematics(second,name)))
                throw std::runtime_error("Unsupported TurbulenceTransport mathematics for "+name
                    +"; the frozen kernel implements source/diffusion with no convection.");
        }
        const auto& closure=system.state.at("mu_t");
        if (closure.role!=StateRole::Derived || closure.storageBinding!=StorageBinding::ProviderDistributed
            || closure.storageKey!="mu_t")
            throw std::runtime_error("TurbulenceTransport requires the provider-owned mu_t closure.");
        CompiledEquationCall result;
        result.physicalResources.push_back({step.occurrence,"model.turbulence",typeid(FDM::ITransportModel),true});
        result.source=step;result.target=bindTarget(system,step);
        result.equationMethod=std::string(id());
        result.backendOperation=OpIds::TurbulenceAdvance;
        result.calls={{step.equation,step.target.symbol}};
        result.writes={step.target.symbol,"mu_t"};
        result.reads={"rho","U","k",second,"mu_t"};
        result.stateUses={{"rho"},{"U"},{"k"},{second}};
        result.coefficientProvides={"turbulence.mu_t.step"};result.coefficientStates={"mu_t"};
        result.coefficientVersion=StateVersion::Frozen;
        result.stateEffects={{"k",false,true},{second,false,true},{"mu_t",true,true},{"turbulence.mu_t.step",true,false}};
        result.capabilityRequirements={{"velocity.immersed","turbulence.boundary.immersed",
            "requires matching immersed mean-flow wall, wall distance, turbulence boundary and velocity time-level treatment",
            "momentum","rhoU","immersed.ghost"}};
        if (second=="omega" && system.immersed
            && system.immersed->wallClosure==FDM::ImmersedWallClosure::StationaryNoSlipAdiabatic) {
            if (!system.state.contains("immersed.wallDistance")
                || !producedBy(system.state.at("immersed.wallDistance"),"geometry.immersedWallDistance"))
                throw std::runtime_error("SST immersed wall requires the original geometry wall-distance view.");
            result.stateUses.push_back({"immersed.wallDistance"});
            result.boundCapabilities={{"turbulence.boundary.immersed","momentum","rhoU","immersed.ghost"}};
            result.capabilityRequirements.push_back({{},"meanflow.boundary.immersed.viscous",
                "SST wall closure must bind the same stationary no-slip mean-flow wall",
                "momentum","rhoU","immersed.ghost"});
        }
        result.fusionKey=std::string(id())+"/"+model;
        result.fusionMembers={{"k","k"},{second,second}};
        result.temporalMethod="physical-step explicit in-place update";
        result.requirements={"single-fluid serial single patch","complete adjacent RAS pair",
            "boundary -> WriteOwned -> ReadHalo -> correction -> boundary -> WriteOwned -> ReadHalo"};
        if (second=="omega" && system.immersed
            && system.immersed->wallClosure==FDM::ImmersedWallClosure::StationaryNoSlipAdiabatic)
            result.requirements.push_back("SST consumes boundary-ready step-entry velocity; mu_t frozen for this physical step");
        result.operations.push_back({OpIds::TurbulenceAdvance,"advance frozen RAS pair once",
            OperationStage::Prepare,{OperationCapability::SingleFluidTurbulenceTransport},
            {OriginKind::Generated,std::string(id())}});
        return result;
    }
};

// Smagorinsky refresh is an algebraic closure, never a transported RAS pair.
class TurbulenceClosureMethod final : public IProvider {
public:
    std::string_view id() const override { return "TurbulenceClosure"; }
    std::string_view runtimeProvider() const override { return "flow.turbulence-closure"; }
    CompiledEquationCall compile(const ExecutableEquationSystem& system,
                                const EquationCall& step,const NumericalBinding&) const override {
        const Equation expected{"mu_t",Turbulence::mathematicalSymbol("mu_t"),
            Turbulence::mathematicalSymbol("Smagorinsky.eddyDynamicViscosity"),{}};
        if (step.equation!="mu_t" || step.target.symbol!="mu_t" || step.target.kind!=TargetKind::Physical
            || canonicalFormula(system.registry.at(step.equation))!=canonicalFormula(expected)
            || system.state.at("mu_t").storageBinding!=StorageBinding::ProviderDistributed)
            throw std::runtime_error("Unsupported TurbulenceClosure equation/state contract.");
        CompiledEquationCall result;
        result.physicalResources.push_back({step.occurrence,"model.turbulence",typeid(FDM::ITransportModel),true});
        result.source=step;result.target=bindTarget(system,step);
        result.equationMethod=std::string(id());
        result.backendOperation=OpIds::TurbulenceClosureRefresh;
        result.calls={{"mu_t","mu_t"}};result.writes={"mu_t"};result.reads={"rho","U"};
        result.coefficientProvides={"turbulence.mu_t.step"};result.coefficientStates={"mu_t"};
        result.coefficientVersion=StateVersion::Frozen;
        result.stateUses={{"rho"},{"U"}};result.stateEffects={{"mu_t",true,true},{"turbulence.mu_t.step",true,false}};
        result.capabilityRequirements={{"velocity.immersed","turbulence.boundary.immersed",
            "immersed strain-rate and wall/boundary closure treatment are unavailable"}};
        result.operations.push_back({OpIds::TurbulenceClosureRefresh,"refresh algebraic eddy viscosity",
            OperationStage::Prepare,{OperationCapability::SingleFluidTurbulenceClosure},{}});
        return result;
    }
};

class DirectEvaluationMethod final : public IProvider {
public:
    std::string_view id() const override { return "DirectEvaluation"; }
    CompiledEquationCall compile(const ExecutableEquationSystem& system,
                                const EquationCall& step,
                                const NumericalBinding&) const override {
        requireOutput(system,step);
        if (step.target.symbol.empty())
            throw std::runtime_error("DirectEvaluation requires one primary Output symbol.");
        const auto& formula=system.registry.at(step.equation);
        const auto& symbol=step.target.symbol;
        if (formula.lhs.kind()!=FormulaExpr::Kind::Symbol
            || formula.lhs.name()!=symbol)
            throw std::runtime_error("DirectEvaluation Equation '"+formula.id
                +"' does not assign the declared primary Output '"+symbol+"'.");
        CompiledEquationCall result;
        result.source=step;
        result.target=bindTarget(system,step);
        result.equationMethod=std::string(id());
        result.backendOperation="formula.direct.evaluate."+step.equation;
        result.calls={{formula.id,symbol}};
        result.writes={symbol};
        if (!result.target.workspace.empty())
            result.workspaceProvides={result.target.workspace};
        collectSymbols(formula.rhs,result.reads);
        return result;
    }
};

class LinearEquationMethod final : public IProvider {
public:
    std::string_view id() const override { return "LinearEquation"; }
    CompiledEquationCall compile(const ExecutableEquationSystem& system,
                                const EquationCall& step,
                                const NumericalBinding&) const override {
        requireOutput(system,step);
        if (step.target.symbol.empty())
            throw std::runtime_error("LinearEquation requires one primary Output symbol.");
        const auto& formula=system.registry.at(step.equation);
        const auto& symbol=step.target.symbol;
        if(inspectFormula(formula.lhs,symbol).targetDegree>1 || inspectFormula(formula.rhs,symbol).targetDegree>1)
            throw std::runtime_error("Unsupported: LinearEquation has no nonlinear/Jacobian implementation for target "+symbol);
        if (!dependsOn(formula.lhs,symbol)
            && !dependsOn(formula.rhs,symbol))
            throw std::runtime_error("LinearEquation Equation '"+formula.id
                +"' does not contain its declared Output '"+symbol+"'.");
        CompiledEquationCall result;
        result.source=step;
        result.target=bindTarget(system,step);
        result.equationMethod=std::string(id());
        result.backendOperation="formula.linear.solve."+step.equation;
        result.calls={{formula.id,symbol}};
        result.writes={symbol};
        if (!result.target.workspace.empty())
            result.workspaceProvides={result.target.workspace};
        result.requirements={"linear assembly provider","GlobalDofSystem",
                             "linear solver backend"};
        collectSymbols(formula.lhs,result.reads);
        collectSymbols(formula.rhs,result.reads);
        return result;
    }
};

/// Pressure methods compile only the local numerical topology. The existing
/// PressureOperators callbacks retain gauge, boundary and canonical-face math.
std::vector<ExecutableOperation> pressureOperations() {
    std::vector<ExecutableOperation> operations;
        const Provenance source{OriginKind::Generated,"pressureConstraint"};
        const auto declare = [&](OperationStage stage, const char* id,
                                 const char* name,
                                 std::vector<OperationCapability> needs) {
            operations.push_back(
                {id,name,stage,std::move(needs),source});
        };
        declare(OperationStage::Prepare,OpIds::PressurePrepare,
                "prepare pressure schedule",{OperationCapability::PressureSchedule});
        declare(OperationStage::FixedTimeStepBegin,OpIds::PressureStepBegin,
                "freeze physical-time base state",
                {OperationCapability::FixedTimeIteration});
        declare(OperationStage::IterationBegin,OpIds::PressureIterationBegin,
                "begin fixed-time outer iterate",
                {OperationCapability::FixedTimeIteration});
        declare(OperationStage::MomentumAssemble,OpIds::MomentumAssemble,
                "assemble momentum predictor",{OperationCapability::MomentumPredictor});
        declare(OperationStage::MomentumSolve,OpIds::MomentumSolve,
                "advance momentum predictor",{OperationCapability::MomentumPredictor});
        declare(OperationStage::PressureBoundaryPrepare,
                OpIds::PressureBoundaryPrepare,"prepare pressure boundary state",
                {OperationCapability::PressureBoundary});
        declare(OperationStage::PressureAssemble,OpIds::PressureAssemble,
                "assemble pressure correction",{OperationCapability::PressureCorrection,
                                                  OperationCapability::PressureLinearSolve});
        declare(OperationStage::PressureSolve,OpIds::PressureSolve,
                "solve pressure correction",{OperationCapability::PressureCorrection,
                                               OperationCapability::PressureLinearSolve});
        declare(OperationStage::PressureUpdatePrepare,
                OpIds::PressureUpdatePrepare,"prepare pressure update",
                {OperationCapability::PressureCorrection});
        declare(OperationStage::VelocityCorrect,OpIds::VelocityCorrect,
                "velocity correction",{OperationCapability::VelocityCorrection});
        declare(OperationStage::FluxCorrect,OpIds::FluxCorrect,
                "refresh derived face-flux state",{OperationCapability::FluxCorrection});
        declare(OperationStage::CorrectionCommit,
                OpIds::PressureCorrectionCommit,
                "commit pressure-corrected state",{OperationCapability::PressureCorrection});
        declare(OperationStage::RelaxationApply,OpIds::PressureRelaxationApply,
                "relax fixed-point solution",
                {OperationCapability::FixedTimeIteration});
        declare(OperationStage::FluxConsistencyRestore,
                OpIds::PressureFluxConsistencyRestore,
                "restore face flux for the relaxed iterate",
                {OperationCapability::FixedTimeIteration,
                 OperationCapability::FluxCorrection});
        declare(OperationStage::ConvergenceEvaluate,
                OpIds::PressureConvergenceEvaluate,
                "evaluate fixed-point residual",
                {OperationCapability::FixedTimeIteration});
        declare(OperationStage::IterationEnd,OpIds::PressureIterationEnd,
                "end fixed-time outer iterate",
                {OperationCapability::FixedTimeIteration});
        declare(OperationStage::StepCommit,OpIds::PressureStepCommit,
                "commit corrected state",{OperationCapability::PressureSchedule});
    return operations;
}

// These resources alias the existing pressure operator arrays. Publication is
// the authored correction tail, not a read-side refresh of physical STATE.
void pressureLifecycleContract(CompiledEquationCall& call,bool conservative) {
    const auto& method=call.equationMethod;
    if (method=="VelocityCorrection" || method=="PressureUpdate" || method=="FluxCorrection"
        || method=="ConservativeVelocityCorrection" || method=="ConservativePressureUpdate" || method=="ConservativeFluxCorrection") {
        CompiledStateView view;view.symbol="deltaP";view.kind=StateViewKind::Correction;call.viewReads.push_back(view);
    }
    if (method=="ConservativePressurePublication") {
        CompiledStateView view;view.symbol="p";view.kind=StateViewKind::Working;call.viewReads.push_back(view);
    }
    call.commitRequires={"pressure.published"};
    const bool predictor=method=="PressureMomentum" || method=="ConservativePressureMomentum";
    if (predictor) {
        call.workspaceInvalidates={"pressureCorrection","preparedPressure","derivedFluxValidity","pressure.published","rAU","faceResponse","correctedFlux"};
        if (conservative) {
            call.stateUses={{"rho"},{"U"},{"p"}};
            call.stateEffects={{"rho"},{"rhoU"},{"rhoE"}};
        } else {
            call.stateUses={{"U"},{"p"}};
            if (call.target.kind==TargetKind::Physical) call.stateEffects={{"U",true,true}};
            call.workspaceProvides.push_back("pressure.working.U");
            call.requirements.push_back("fixed-time predictor reads original baseVelocity_ (OldTime); iteration velocity is Working(U)");
        }
    } else if (method=="PressureCorrection") {
        call.stateUses={{"faceResponse",StateVersion::Frozen,false,"faceResponse"},
            {"correctedFlux",StateVersion::Frozen,false,"correctedFlux"}};
        call.workspaceRequires.push_back("pressure.working.U");
    } else if (method=="ConservativeVelocityCorrection") {
        call.stateEffects={{"rhoU"}};
        call.workspaceInvalidates={"pressure.published"};
    } else if (method=="VelocityCorrection") {
        call.workspaceRequires.push_back("pressure.working.U");
        call.stateUses={{"rAU",StateVersion::Frozen,false,"rAU"}};
        call.workspaceInvalidates={"pressure.published"};
        call.stateEffects={{"U",true,true}};
    } else if (method=="PressureUpdate" || method=="ConservativePressureUpdate") {
        call.workspaceInvalidates={"pressure.published"};
        if (!conservative) call.stateEffects={{"p",true,true}};
    } else if (method=="FluxCorrection" || method=="ConservativePressurePublication") {
        call.workspaceInvalidates={"pressureCorrection"};
        if (!conservative) {
            call.workspaceInvalidates.push_back("correctedFlux");
            call.workspaceProvides.push_back("correctedFlux");
        }
        call.workspaceProvides.push_back("pressure.published");
        call.stepWorkspaces.push_back("pressure.published");
        if (!conservative) {
            call.stateEffects={{"phi",true,true}};
        } else call.stateEffects={{"rhoE",true,false}};
    } else if (method=="FixedTimeRelaxation" || method=="ConservativeFixedTimeRelaxation") {
        call.workspaceRequires.push_back("pressure.published");
        call.workspaceInvalidates={"pressure.published"};
        if (!conservative) call.stateEffects={{"U",true,true},{"p",true,true}};
    } else if (method=="FluxConsistency" || method=="ConservativeFluxConsistency") {
        call.workspaceProvides.push_back("pressure.published");call.stepWorkspaces.push_back("pressure.published");
        if (!conservative) call.stateEffects={{"phi",true,true}};
    } else if (method=="ResidualConvergence" || method=="ConservativeResidualConvergence") {
        call.workspaceRequires.push_back("pressure.published");
    }
}

// This implementation consumes additive primitive momentum terms only. Signs,
// repeated stencil terms and non-symbol coefficients need different kernels.
void lowerPressureMomentum(const ExecutableEquationSystem& system,const CompiledEquationCall& call,
        const FDM::NumericalRecipeSet& recipes,const TermProviderCatalog& catalog,
        const std::vector<FormulaOperatorBinding>& bindings,CompiledNumericalSystem& result) {
    const auto& equation=system.registry.at(call.source.equation);
    const auto context="PressureMomentum step="+call.source.occurrence+" equation="+equation.id+" target="+call.target.symbol;
    const auto fail=[&](const std::string& path,const std::string& reason) {
        throw std::runtime_error(context+" occurrence="+path+": "+reason);
    };
    if(result.pressureMomentum)fail("HOW","multiple predictors require independently bound pressure workspaces");
    CompiledPressureMomentum momentum{equation.id,call.source.occurrence,{},{}};
    int transient=0,convection=0,diffusion=0;
    const auto visit=[&](const auto& self,const FormulaExpr& e,CompiledSpatialBinding::Side side,const std::string& path)->void {
        using K=FormulaExpr::Kind;using Side=CompiledSpatialBinding::Side;
        const auto occurrence=e.occurrence().empty()?path:e.occurrence();
        if(e.kind()==K::Add) {self(self,e.arguments()[0],side,path+"/0");self(self,e.arguments()[1],side,path+"/1");return;}
        if(e.kind()==K::Constant && e.constant()==0.0 && side==Side::Right)return;
        if(e.kind()!=K::Operator)fail(occurrence,"requires positive additive ddt(U), div(momentumFlux), diffusion(nu,U), or RHS source(symbol); no scaling/sign lowering");
        const auto symbol=[&](std::size_t n,std::string_view expected) {
            return e.arguments().size()>n && e.arguments()[n].kind()==K::Symbol && e.arguments()[n].name()==expected;
        };
        const FDM::TermRecipe* requested=nullptr;
        std::string primary,secondary;
        if(e.operation()==Math::Op::Ddt) {
            if(side!=Side::Left || e.arguments().size()!=1 || !symbol(0,"U") || ++transient!=1)
                fail(occurrence,"requires exactly one unscaled LHS ddt(U)");
            if(selectFormulaBinding(equation,e,occurrence,bindings))
                fail(occurrence,"ddt uses the unique selected time recipe; no spatial operator override is implemented");
            momentum.timeOccurrence=occurrence;return;
        } else if(e.operation()==Math::Op::Div) {
            if(side!=Side::Left || e.arguments().size()!=1 || !symbol(0,"momentumFlux") || ++convection!=1)
                fail(occurrence,"requires one LHS div(momentumFlux); repeated convection is not accumulated by this predictor");
            if(!recipes.convection || recipes.convection->id()!=FDM::TermRecipeId::PrimitiveUpwind1)
                fail(occurrence,"selected pressure stencil requires PrimitiveUpwind1");
            requested=&*recipes.convection;primary="momentumFlux";
        } else if(e.operation()==Math::Op::Diffusion) {
            if(side!=Side::Left || e.arguments().size()!=2 || !symbol(0,"nu") || !symbol(1,"U") || ++diffusion!=1)
                fail(occurrence,"requires at most one LHS diffusion(nu,U); variable/scaled coefficients are not implemented");
            if(!recipes.diffusion || recipes.diffusion->id()!=FDM::TermRecipeId::Central2Explicit)
                fail(occurrence,"selected pressure stencil requires Central2Explicit; no Central4 implementation is bound");
            if(!system.state.contains("nu") || !producedBy(system.state.at("nu"),"closure.thermodynamic.nu"))
                fail(occurrence,"nu must use the existing constant transport projection; this predictor has no custom coefficient field/producer port");
            requested=&*recipes.diffusion;primary="U";secondary="nu";
        } else if(e.operation()==Math::Op::Source) {
            if(side!=Side::Right || e.arguments().size()!=1)
                fail(occurrence,"primitive source requires one RHS source expression returning acceleration");
            primary=e.arguments()[0].kind()==K::Symbol?e.arguments()[0].name():"expression";
        } else fail(occurrence,std::string("unsupported pressure operator ")+Math::name(e.operation()));
        const auto selected=catalog.resolve({equation,e,occurrence,"U",call.equationMethod,requested,&system.state,&call.stateUses},bindings);
        if(selected.status!=BindingStatus::Resolved)fail(occurrence,selected.reason);
        if(requested && (!selected.recipe || selected.recipe->id()!=requested->id() || !selected.primitiveSpatial))
            fail(occurrence,"selected spatial implementation has no compatible primitive momentum kernel: "+selected.id);
        if(!requested && !selected.primitiveSource && !selected.fieldSource.bind)
            fail(occurrence,"selected source has no primitive acceleration callback: "+selected.id);
        momentum.terms.push_back(result.operators.size());
        result.operators.emplace_back(equation.id,call.source.equation,"U",call.equationMethod,occurrence,e.operation(),
            primary,secondary,side,selected.recipe,selected.id,selected.owner,selected.compiledDataAvailable,
            selected.primitiveSource,ConservativeSourceKernel{},selected.primitiveSpatial);
        result.operators.back().fieldSource=selected.fieldSource;
        if(std::find(result.providerRequirements.begin(),result.providerRequirements.end(),selected.id)==result.providerRequirements.end())
            result.providerRequirements.push_back(selected.id);
        if(selected.recipe)result.recipeBindings.push_back({*selected.recipe,RecipeConsumerKind::EquationTerm,equation.id+"@"+occurrence});
        if(requested)result.requiredHaloWidth=std::max(result.requiredHaloWidth,requested->haloWidth());
    };
    visit(visit,equation.lhs,CompiledSpatialBinding::Side::Left,"lhs");
    visit(visit,equation.rhs,CompiledSpatialBinding::Side::Right,"rhs");
    if(transient!=1 || convection!=1)fail("equation","requires one ddt(U) and one div(momentumFlux)");
    if(std::find(result.workspaceRequirements.begin(),result.workspaceRequirements.end(),call.residualWorkspace)==result.workspaceRequirements.end())
        result.workspaceRequirements.push_back(call.residualWorkspace);
    result.pressureMomentum=std::move(momentum);
}

// The correction kernels still implement these finite structures, not a
// general AST evaluator. Reject every extra term until its kernel is replaced.
void requirePressureRelation(const Equation& equation,std::string_view method) {
    using K=FormulaExpr::Kind;
    const auto shape=[&](const auto& self,const FormulaExpr& e)->std::string {
        if(e.kind()==K::Symbol)return e.name();
        if(e.kind()==K::Operator && e.operation()==Math::Op::Extension
            && e.extensionRule()!=PressureMathematics::op(e.name(),e.arguments()).extensionRule())
            throw std::runtime_error("Pressure kernel requires its declared extension semantics: "+e.name());
        std::string name;
        switch(e.kind()) {
        case K::Add:name="add";break;case K::Subtract:name="sub";break;
        case K::Multiply:name="mul";break;case K::Divide:name="divide";break;
        case K::Negate:name="neg";break;case K::Operator:name=e.name();break;
        default:return "literal";
        }
        name+="(";bool first=true;
        for(const auto& a:e.arguments()) {if(!first)name+=",";first=false;name+=self(self,a);}
        return name+")";
    };
    std::string_view expected;
    if(method=="PressureCorrection")expected="neg(div(mul(faceResponse,grad(deltaP))))=neg(faceVolumeDivergence(correctedFlux))";
    else if(method=="PressureUpdate")expected="p=add(p,deltaP)";
    else if(method=="VelocityCorrection")expected="U=sub(U,mul(rAU,grad(deltaP)))";
    else if(method=="FluxCorrection")expected="phi=sub(phi,pressureFlux(faceResponse,deltaP))";
    else if(method=="FixedTimeRelaxation")expected="iterate=relax(iterate,laggedIterate)";
    else if(method=="FluxConsistency")expected="phi=consistentFlux(U,p)";
    else if(method=="ResidualConvergence")expected="converged=residualConvergence(iterate,laggedIterate)";
    else return;
    if(shape(shape,equation.lhs)+"="+shape(shape,equation.rhs)!=expected)
        throw std::runtime_error("Unsupported fixed pressure relation: equation="+equation.id+" method="+std::string(method)+"; kernel cannot consume the changed AST");
}

class PressureMethod final {
public:
    enum class OutputContract { Assignment, Relation };
    PressureMethod(std::string_view method,std::string_view output,
                   OutputContract contract,std::vector<SolvePlanNode> leaves,
                   std::vector<std::string> needs,
                   std::vector<std::string> provides,
                   std::vector<std::string> requirements,
                   std::size_t sourceInputs=0)
        : id_(method),output_(output),contract_(contract),leaves_(std::move(leaves)),
          needs_(std::move(needs)),provides_(std::move(provides)),
          requirements_(std::move(requirements)),sourceInputs_(sourceInputs) {}
    std::string_view id() const { return id_; }
    std::string_view runtimeProvider() const {
        return id_=="FluxCorrection" ? "flow.rhie-chow" : "flow.pressure-operators";
    }
    std::string_view operationProvider(const OpId& operation) const {
        // Flux publication uses the pressure backend's existing state authority.
        return operation==OpIds::FluxCorrect ? "flow.rhie-chow" : "flow.pressure-operators";
    }
    CompiledEquationCall compile(const ExecutableEquationSystem& system,
                                const EquationCall& step,
                                const NumericalBinding& binding) const {
        requireOutput(system,step);
        const bool supportedTarget=id_=="PressureMomentum"
            ? (step.target.kind==TargetKind::Physical || step.target.kind==TargetKind::Working)
            : step.target.kind==(id_=="PressureCorrection" ? TargetKind::Correction
                : TargetKind::Physical);
        if (!supportedTarget)
            throw std::runtime_error("Unsupported: provider "+std::string(id_)+" cannot realize semantic target "+targetText(step.target));
        if (step.target.symbol.empty() || step.target.symbol!=
            output_)
            throw std::runtime_error("EquationMethod '"+std::string(id_)
                +"' requires explicit Output '"+std::string(output_)+"'.");
        const auto& formula=system.registry.at(step.equation);
        requirePressureRelation(formula,id_);
        if(!binding.parameters.empty())throw std::runtime_error("Pressure method has no per-call scalar parameter implementation: "+step.occurrence);
        if (contract_==OutputContract::Relation) {
            if (!dependsOn(formula.lhs,output_)
                && !dependsOn(formula.rhs,output_))
                throw std::runtime_error("Pressure method Equation does not contain its Output.");
        } else if (formula.lhs.kind()!=FormulaExpr::Kind::Symbol
                   || formula.lhs.name()!=output_) {
            throw std::runtime_error("Pressure method Equation does not assign its Output.");
        }
        if (binding.inputs.size()!=sourceInputs_)
            throw std::runtime_error("EquationMethod '"+std::string(id_)
                +"' has missing or unexpected source-math inputs.");
        for (const auto& input:binding.inputs) {
            if (std::none_of(system.registry.entries().begin(),
                    system.registry.entries().end(),
                    [&](const Equation& candidate) { return candidate.id==input; }))
                throw std::runtime_error("EquationMethod source mathematics '"
                    +input+"' is unavailable.");
        }
        CompiledEquationCall result;
        result.source=step;result.backendProvider=runtimeProvider();
        result.target=bindTarget(system,step);
        result.equationMethod=std::string(id_);
        result.target=bindTarget(system,step,provides_.empty() ? "" : provides_.back());
        result.target.viewOwner=StateViewOwner::NumericalProvider;
        result.temporalCapabilities={FDM::TimeRecipeId::ForwardEuler};
        if (id_=="PressureMomentum") {
            const auto& density=system.state.at("rho");
            if(!density.constantValue || !std::isfinite(*density.constantValue) || *density.constantValue<=0.0)
                throw std::runtime_error("PressureMomentum requires positive constant STATE rho: step="+step.occurrence+" equation="+step.equation);
            for(const auto& read:step.inputs)
                if(read.version!=StateVersion::Current
                    || read.kind!=TargetKind::Physical || read.halo || !read.snapshot.empty() || read.perIteration)
                    throw std::runtime_error("PressureMomentum unsupported input version/port: step="+step.occurrence+" equation="+step.equation+" input="+read.symbol+"; requires current predictor U/p, original fixed-time base is internal workspace");
            result.stepInputs={{OpIds::PressurePrepare,std::string(runtimeProvider())}};
            result.physicalPublication={OpIds::PressureStepCommit,std::string(runtimeProvider())};
            if(step.target.kind==TargetKind::Working) {
                result.stepInputs.push_back({OpIds::PressureStepBegin,std::string(runtimeProvider())});
                result.iterationInputs={{OpIds::PressureIterationBegin,std::string(runtimeProvider())}};
                result.iterationOutputs={{OpIds::PressureIterationEnd,std::string(runtimeProvider())}};
                result.iterationSignalRequired=true;
            }
            result.temporalMethod="fixed physical dt; U_n + dt R(U_iterate), pressure constraint correction";
            result.spatialTerms=true;
            result.primitiveSourceRequired=true;
            result.residualWorkspace="pressureMomentumWorkspace";
            result.operations=pressureOperations();
        }
        result.calls={{formula.id,std::string(output_)}};
        result.writes={std::string(output_)};
        result.requirements=requirements_;
        result.workspaceRequires=needs_;
        result.workspaceProvides=provides_;
        result.stepWorkspaces=provides_;
        result.sourceMathInputs=binding.inputs;
        result.reads=binding.inputs;
        collectSymbols(formula.lhs,result.reads);
        collectSymbols(formula.rhs,result.reads);
        result.fragment.kind=PlanNodeKind::Sequence;
        result.fragment.id=std::string(id_)+".fragment";
        result.fragment.name=std::string(id_)+" numerical fragment";
        pressureLifecycleContract(result,false);
        if(id_=="PressureMomentum")bindSourceInputs(system,step,formula,result);
        result.fragment.children=leaves_;
        for(auto& leaf:result.fragment.children)leaf.provider=operationProvider(leaf.operation);
        result.fragment.children.back().equationCalls=result.calls;
        return result;
    }
private:
    std::string_view id_;
    std::string_view output_;
    OutputContract contract_;
    std::vector<SolvePlanNode> leaves_;
    std::vector<std::string> needs_;
    std::vector<std::string> provides_;
    std::vector<std::string> requirements_;
    std::size_t sourceInputs_;
};

// The conservative predictor publishes the packed physical state in the
// frozen algorithm. Its fusion contract explicitly includes mass and energy;
// Working(p) and Correction(p) alias the Corrector's existing workspace.
class ConservativePressureMethod final : public IProvider {
public:
    void lowerSpatial(const ExecutableEquationSystem& system,const CompiledEquationCall& call,
            const FDM::NumericalRecipeSet& recipes,const TermProviderCatalog& catalog,
            const std::vector<FormulaOperatorBinding>& bindings,CompiledNumericalSystem& result) const override {
        FusedSpatial::lower(system,call,recipes,catalog,bindings,result);
    }
    bool usesSpatialRecipes() const override { return true; }
    ConservativePressureMethod(std::string_view method,std::string_view equation,
            std::string_view output,TargetKind kind,std::vector<const char*> operations,
            std::vector<std::string> needs={},std::vector<std::string> provides={})
        : method_(method),equation_(equation),output_(output),kind_(kind),
          operations_(std::move(operations)),needs_(std::move(needs)),provides_(std::move(provides)) {}
    std::string_view id() const override { return method_; }
    std::string_view runtimeProvider() const override { return "flow.conservative"; }
    CompiledEquationCall compile(const ExecutableEquationSystem& system,
            const EquationCall& call,const NumericalBinding& binding) const override {
        requireOutput(system,call);
        if (call.target.symbol!=output_ || call.target.kind!=kind_)
            throw std::runtime_error("Unsupported: provider "+std::string(method_)+" requires "
                +std::string(equation_)+" -> "+std::string(output_)+" with its declared semantic target.");
        CompiledEquationCall result;
        result.source=call;result.equationMethod=method_;
        result.target=bindTarget(system,call,provides_.empty() ? "" : provides_.back());
        result.target.viewOwner=StateViewOwner::NumericalProvider;
        result.temporalCapabilities={FDM::TimeRecipeId::ForwardEuler};
        result.workspaceRequires=needs_;result.workspaceProvides=provides_;
        result.stepWorkspaces=provides_;
        const auto& formula=system.registry.at(call.equation);
        if (method_=="ConservativePressureMomentum") {
            result.stepInputs={{OpIds::PressurePrepare,std::string(runtimeProvider())}};
            result.physicalPublication={OpIds::PressureStepCommit,std::string(runtimeProvider())};
            result.temporalMethod="one forwardEuler conservative predictor; fixed-time pressure corrections";
            if (binding.inputs.size()!=2)
                throw std::runtime_error("ConservativePressureMomentum requires explicit continuity/energy fusion inputs.");
            result.spatialTerms=true;result.residualWorkspace="residual";
            result.operations=pressureOperations();
            result.sourceMathInputs=binding.inputs;
            for (const auto& math:std::vector<CompiledMathRef>{{binding.inputs[0],"rho"},{call.equation,"rhoU"},{binding.inputs[1],"rhoE"}}) {
                const auto& input=system.registry.at(math.equation);
                if (!containsDdt(input.lhs,math.target))
                    throw std::runtime_error("Conservative pressure predictor requires ddt("+math.target+").");
                result.calls.push_back(math);result.writes.push_back(math.target);
                collectSymbols(input.lhs,result.reads);collectSymbols(input.rhs,result.reads);
            }
            result.requirements={"five-component packed conservative state","EOS closure",
                "boundary-ready state","single-stage fused conservative residual"};
        } else {
            if (!binding.inputs.empty()) throw std::runtime_error("Unexpected conservative pressure source-math inputs.");
            const auto relations=conservativePressureRelations();
            const auto expected=std::find_if(relations.begin(),relations.end(),[&](const auto& value) { return value.id==equation_; });
            if (expected==relations.end() || canonicalFormula(formula)!=canonicalFormula(*expected))
                throw std::runtime_error("Unsupported: selected conservative pressure provider cannot execute changed mathematics: "+call.equation);
            result.calls={{formula.id,std::string(output_)}};result.writes={std::string(output_)};
            collectSymbols(formula.lhs,result.reads);collectSymbols(formula.rhs,result.reads);
            result.requirements={"bound conservative pressure relation","EOS closure","pressure gauge/reference"};
        }
        result.fragment.kind=PlanNodeKind::Sequence;
        if(method_=="ConservativeFixedTimeRelaxation") {
            result.iterationSignalRequired=true;
            result.stepInputs={{OpIds::PressureStepBegin,std::string(runtimeProvider())}};
            result.iterationInputs={{OpIds::PressureIterationBegin,std::string(runtimeProvider())}};
            result.iterationOutputs={{OpIds::PressureIterationEnd,std::string(runtimeProvider())}};
        }
        result.fragment.id=std::string(method_)+".fragment";
        result.fragment.name=std::string(method_)+" numerical fragment";
        pressureLifecycleContract(result,true);
        for (const auto* operation:operations_)
            result.fragment.children.push_back(methodLeaf(PlanNodeKind::Update,
                std::string(method_)+"."+operation,operation,operation,runtimeProvider()));
        result.fragment.children.back().equationCalls=result.calls;
        return result;
    }
private:
    std::string_view method_,equation_,output_;
    TargetKind kind_;
    std::vector<const char*> operations_;
    std::vector<std::string> needs_,provides_;
};

class PressureRelation final {
public:
    PressureRelation(std::string id,std::string operation)
        : id_(std::move(id)),operation_(std::move(operation)) {}
    std::string_view id() const { return id_; }
    std::string_view runtimeProvider() const { return "flow.pressure-operators"; }
    CompiledEquationCall compile(const ExecutableEquationSystem& system,
            const EquationCall& call,const NumericalBinding& binding) const {
        requireOutput(system,call);
        if(!binding.inputs.empty() || !binding.parameters.empty() || !binding.operators.empty())
            throw std::runtime_error("Fixed pressure relation has no auxiliary input/parameter/operator implementation: "+call.occurrence);
        const auto& equation=system.registry.at(call.equation);
        requirePressureRelation(equation,id_);
        CompiledEquationCall result;
        result.source=call; result.equationMethod=id_;result.backendProvider=runtimeProvider();
        result.target=bindTarget(system,call);
        result.target.viewOwner=StateViewOwner::NumericalProvider;
        result.temporalCapabilities={FDM::TimeRecipeId::ForwardEuler};
        result.backendOperation=operation_;
        pressureLifecycleContract(result,false);
        result.calls.push_back({equation.id,call.target.symbol});
        result.writes.push_back(call.target.symbol);
        if (!result.target.workspace.empty()) result.workspaceProvides.push_back(result.target.workspace);
        collectSymbols(equation.lhs,result.reads);
        collectSymbols(equation.rhs,result.reads);
        return result;
    }
private:
    std::string id_,operation_;
};

} // namespace

TemporalMethodRegistry builtinTemporalMethods() {
    static const ForwardEulerMethod euler;
    static const SSPRK3Method ssp;
    static const ClassicalRK4Method rk4;
    TemporalMethodRegistry result;
    result.add(euler);
    result.add(ssp);
    result.add(rk4);
    return result;
}

std::optional<CompiledEquationCall> compileDirectEquation(const ExecutableEquationSystem& system,
        const EquationCall& call,const NumericalBinding& binding) {
    if(binding.method=="LocalAffine")return compileLocalAffine(system,call,binding);
    if(binding.method=="ConservativeResidual") {
        auto result=compileConservativeResidual(system,call,binding);result.backendProvider="flow.conservative";return result;
    }
    static const PressureMethod momentum("PressureMomentum","U",PressureMethod::OutputContract::Relation,
        {methodLeaf(PlanNodeKind::Assemble,"PressureMomentum.assemble",
            "assemble momentum predictor",OpIds::MomentumAssemble),
         methodLeaf(PlanNodeKind::Solve,"PressureMomentum.solve",
            "solve momentum predictor",OpIds::MomentumSolve)},
        {},{"HbyA","rAU","faceResponse","correctedFlux","pressureMomentumWorkspace"},
        {"bound Equation occurrence mathematics",
         "velocity boundary","RhieChow face prediction","canonical face COPY"});
    static const PressureMethod correction("PressureCorrection","deltaP",
        PressureMethod::OutputContract::Relation,
        {methodLeaf(PlanNodeKind::Update,"PressureCorrection.boundary",
            "prepare pressure boundary state",OpIds::PressureBoundaryPrepare),
         methodLeaf(PlanNodeKind::Assemble,"PressureCorrection.assemble",
            "assemble pressure correction",OpIds::PressureAssemble),
         methodLeaf(PlanNodeKind::Solve,"PressureCorrection.solve",
            "solve pressure correction",OpIds::PressureSolve)},
        {"faceResponse","correctedFlux"},{"pressureCorrection"},
        {"pressure boundary","pressure gauge/reference","GlobalDofSystem",
         "linear solver backend"});
    static const PressureMethod pressureUpdate("PressureUpdate","p",
        PressureMethod::OutputContract::Assignment,
        {methodLeaf(PlanNodeKind::Correct,"PressureUpdate.prepare",
            "prepare pressure update",OpIds::PressureUpdatePrepare)},
        {"pressureCorrection"},{},
        {"pressure gauge/reference","pressure boundary"});
    static const PressureMethod velocityCorrection("VelocityCorrection","U",
        PressureMethod::OutputContract::Assignment,
        {methodLeaf(PlanNodeKind::Correct,"VelocityCorrection.correct",
            "correct velocity",OpIds::VelocityCorrect)},
        {"pressureCorrection","rAU"},{},
        {"velocity boundary closure"});
    static const PressureMethod fluxCorrection("FluxCorrection","phi",
        PressureMethod::OutputContract::Assignment,
        {methodLeaf(PlanNodeKind::Correct,"FluxCorrection.correct",
            "correct face flux",OpIds::FluxCorrect),
         methodLeaf(PlanNodeKind::Update,"FluxCorrection.publish",
            "publish corrected pressure state",OpIds::PressureCorrectionCommit)},
        {"pressureCorrection","faceResponse"},{"pressureFaceFlux"},
        {"RhieChow pressureFlux provider","canonical face COPY",
         "continuity defect"});
    static const PressureRelation relax("FixedTimeRelaxation",OpIds::PressureRelaxationApply);
    static const PressureRelation restore("FluxConsistency",OpIds::PressureFluxConsistencyRestore);
    static const PressureRelation convergence("ResidualConvergence",OpIds::PressureConvergenceEvaluate);
    for(const auto* method:{&momentum,&correction,&pressureUpdate,&velocityCorrection,&fluxCorrection})
        if(binding.method==method->id())return method->compile(system,call,binding);
    for(const auto* method:{&relax,&restore,&convergence})
        if(binding.method==method->id())return method->compile(system,call,binding);
    return std::nullopt;
}
bool lowerDirectSpatial(const ExecutableEquationSystem& system,const CompiledEquationCall& call,
        const FDM::NumericalRecipeSet& recipes,const TermProviderCatalog& catalog,
        const std::vector<FormulaOperatorBinding>& bindings,CompiledNumericalSystem& result) {
    if(call.equationMethod=="PressureMomentum")lowerPressureMomentum(system,call,recipes,catalog,bindings,result);
    else if(call.equationMethod=="ConservativeResidual")FusedSpatial::lower(system,call,recipes,catalog,bindings,result);
    else return false;
    return true;
}

std::string validateCompiledConservativeStage(const ExecutableEquationSystem& equations,
    const CompiledSolvePlan& plan) {
    const auto findStage=[&](const auto& self,
                             const System::SolvePlanNode& node)
        -> const System::SolvePlanNode* {
        if (node.operation==System::OpIds::ExplicitStageExecute || node.operation=="flow.stage.rhs") return &node;
        for (const auto& child:node.children)
            if (const auto* found=self(self,child)) return found;
        return nullptr;
    };
    const auto* stage=findStage(findStage,plan.root);
    if (!stage || stage->equationCalls.empty())
        return "Fused explicit stage has no compiled EquationCalls.";
    const auto& compiled=plan.compiledProgram;
    if (!compiled.hasTemporalRoot)
        return "Fused conservative backend requires a TemporalMethod fragment.";
    std::vector<System::CompiledMathRef> calls;
    for (const auto& occurrence:compiled.steps) {
        if (!occurrence.temporalResidual || (!compiled.temporalParticipants.empty() && occurrence.backendProvider!="flow.conservative")) continue;
        if (occurrence.backendOperation!=System::OpIds::ExplicitStageExecute
            || occurrence.backendProvider!="flow.conservative")
            return "Fused explicit stage has an incompatible selected backend.";
        calls.insert(calls.end(),occurrence.calls.begin(),occurrence.calls.end());
    }
    const std::vector<std::string> layout{"rho","rhoU","rhoE"};
    if (calls.size()<layout.size() || stage->equationCalls.size()!=calls.size())
        return "Fused conservative stage requires rho/rhoU/rhoE targets.";
    for (size_t n=0;n<calls.size();++n)
        if ((n<layout.size() && calls[n].target!=layout[n])
            || (n>=layout.size() && std::find(layout.begin(),layout.end(),calls[n].target)==layout.end()
                && (equations.state.at(calls[n].target).storageBinding!=System::StorageBinding::NamedDistributed
                || equations.state.at(calls[n].target).components!=1))
            || stage->equationCalls[n].equation!=calls[n].equation
            || stage->equationCalls[n].target!=calls[n].target)
            return "Fused stage differs from compiled equation occurrence binding.";
    return {};
}


void validateDirectProgram(const ExecutableEquationSystem& equations,const CompiledNumericalSystem& numerics,const CompiledSolvePlan& plan) {
    if(!numerics.pressureMomentum && plan.compiledProgram.temporalParticipants.empty()
        && std::any_of(numerics.operators.begin(),numerics.operators.end(),[](const auto& term){return bool(term.fieldSource.bind);}))
        throw std::runtime_error("Field source requires the implemented synchronous stage binding; the unsplit fused executor has no bound field-reader source port");
    if(numerics.pressureMomentum) {
        const auto reason=validatePressureProgram(equations,numerics,plan);
        if(!reason.empty())throw std::runtime_error(reason);
    }
    if(std::any_of(plan.compiledProgram.steps.begin(),plan.compiledProgram.steps.end(),[](const auto& call){return call.equationMethod=="ConservativeResidual";})) {
        const auto stageContract=validateCompiledConservativeStage(equations,plan);
        if (!stageContract.empty()) throw std::runtime_error(stageContract);
        const std::vector<std::string> layout{"rho","rhoU","rhoE"};
        std::vector<CompiledEquationCall> calls;
        for (const auto& call:plan.compiledProgram.steps) if (call.equationMethod=="ConservativeResidual") calls.push_back(call);
        if (calls.size()!=layout.size()) throw std::runtime_error("ConservativeResidual requires exactly rho/rhoU/rhoE occurrences");
        for (std::size_t i=0;i<layout.size();++i) {
            if (!calls[i].temporalResidual || calls[i].backendProvider!="flow.conservative"
                || calls[i].source.target.symbol!=layout[i])
                throw std::runtime_error("ConservativeResidual has incompatible occurrence/target layout");
            const auto& symbol=equations.state.at(layout[i]);
            const int components=i==1 ? 3 : 1;
            const int offset=i==0 ? 0 : i==1 ? 1 : 4;
            if (symbol.storageBinding!=StorageBinding::PackedDistributed || symbol.storageKey!="conservative"
                || symbol.components!=components || symbol.componentOffset!=offset)
                throw std::runtime_error("ConservativeResidual has incompatible physical storage: "+layout[i]);
        }
    }
}

ProviderRegistry builtinProviders() {
    static const DirectEvaluationMethod evaluation;
    static const LinearEquationMethod linear;
    static const ConservativePressureMethod conservativeMomentum("ConservativePressureMomentum","momentum","rhoU",TargetKind::Physical,
        {OpIds::MomentumAssemble,OpIds::MomentumSolve},{},{"conservativePredictor"});
    static const ConservativePressureMethod conservativeCorrection("ConservativePressureCorrection","pSimple","deltaP",TargetKind::Correction,
        {OpIds::PressureBoundaryPrepare,OpIds::PressureAssemble,OpIds::PressureSolve},{"conservativePredictor"},{"pressureCorrection"});
    static const ConservativePressureMethod conservativeVelocity("ConservativeVelocityCorrection","correctU","rhoU",TargetKind::Physical,
        {OpIds::VelocityCorrect},{"pressureCorrection"});
    static const ConservativePressureMethod conservativeUpdate("ConservativePressureUpdate","correctP","p",TargetKind::Working,
        {OpIds::PressureUpdatePrepare},{"pressureCorrection"},{"preparedPressure"});
    static const ConservativePressureMethod conservativeFlux("ConservativeFluxCorrection","correctFluxp","fluxValidity",TargetKind::Workspace,
        {OpIds::FluxCorrect},{"pressureCorrection"},{"derivedFluxValidity"});
    static const ConservativePressureMethod conservativePublication("ConservativePressurePublication","publishPressure","rhoE",TargetKind::Physical,
        {OpIds::PressureCorrectionCommit},{"preparedPressure","derivedFluxValidity"});
    static const ConservativePressureMethod conservativeRelax("ConservativeFixedTimeRelaxation","relaxIterate","iterate",TargetKind::Workspace,
        {OpIds::PressureRelaxationApply});
    static const ConservativePressureMethod conservativeRestore("ConservativeFluxConsistency","restoreFlux","fluxValidity",TargetKind::Workspace,
        {OpIds::PressureFluxConsistencyRestore},{},{"derivedFluxValidity"});
    static const ConservativePressureMethod conservativeConvergence("ConservativeResidualConvergence","checkConvergence","converged",TargetKind::Workspace,
        {OpIds::PressureConvergenceEvaluate});
    static const TurbulenceTransportMethod turbulence;
    ProviderRegistry result;
    addScalarMethod(result);
    static const SstClosureMethod closure;result.add(closure);
    addEulerianMethods(result);
    addImmersedMethods(result);
    addLevelSetMethods(result);
    addMultiphaseMethods(result);
    addThermalMethod(result);
    result.add(turbulence);
    static const TurbulenceClosureMethod turbulenceClosure;
    result.add(turbulenceClosure);
    for (const auto* provider:{&conservativeMomentum,&conservativeCorrection,&conservativeVelocity,
            &conservativeUpdate,&conservativeFlux,&conservativePublication,&conservativeRelax,
            &conservativeRestore,&conservativeConvergence}) result.add(*provider);
    result.add(evaluation);
    result.add(linear);
    return result;
}


} // namespace SF::System
