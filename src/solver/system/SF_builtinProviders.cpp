#include "core/interfaces/SF_transportModel.h"
#include "SF_fusedSpatial.h"
#include "models/turbulence/SF_transportMathematics.h"
/// Selected numerical implementations compile mathematics; HOW owns the process.
#include "SF_methodObjects.h"
#include "SF_immersedMethods.h"
#include "SF_levelSetMethods.h"
#include "SF_multiphaseMethods.h"
#include "SF_eulerianCoupling.h"
#include "SF_pressureCoupling.h"
#include "SF_affine.h"
#include "SF_block.h"
#include "SF_cell.h"
#include "SF_executionContract.h"
#include "SF_scalarMethod.h"
#include "SF_sstMethod.h"
#include "SF_formulaCompiler.h"
#include "core/system/SF_operationIds.h"
#include <algorithm>
#include <set>
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
                               const std::vector<CompiledMathRef>& calls,
                               const std::vector<TemporalParticipant>& participants,
                               const std::vector<SolvePlanNode>& prefix) {
    if (!body.temporalResidual || calls.empty()
        || body.backendOperation.empty())
        throw std::runtime_error("TemporalMethod requires a compiled residual body and backend.");
    SolvePlanNode root;
    root.kind=PlanNodeKind::Sequence;
    root.id="Explicit.step";
    root.name="explicit time step";
    if (!participants.empty()) {
        const auto leaf=[&](const std::string& op,const std::string& owner,
                const TemporalParticipant* participant=nullptr) {
            auto node=methodLeaf(PlanNodeKind::Update,op,op,op,owner);
            if (participant) {node.occurrence=participant->identity;node.equationCalls=participant->calls;}
            return node;
        };
        for (const auto& p:participants) for (const auto& op:p.preparation)
            root.children.push_back(leaf(op,p.provider,&p));
        root.children.push_back(leaf(TemporalOps::Dt,TemporalOps::Provider));
        root.children.insert(root.children.end(),prefix.begin(),prefix.end());
        for (const auto& p:participants) root.children.push_back(leaf(p.snapshot,p.provider,&p));
        auto stages=methodLeaf(PlanNodeKind::StageLoop,"Explicit.stages","synchronous explicit stages",{});
        stages.repetitions=recipe.stageCount();
        stages.children.push_back(leaf(TemporalOps::Open,TemporalOps::Provider));
        for (const auto& p:participants) stages.children.push_back(leaf(p.prepareStage,p.provider,&p));
        stages.children.push_back(leaf(TemporalOps::Ready,TemporalOps::Provider));
        for(const auto& p:participants)if(!p.coefficientPrepare.empty())stages.children.push_back(leaf(p.coefficientPrepare,p.provider,&p));
        for (const auto& p:participants) stages.children.push_back(leaf(p.rhs,p.provider,&p));
        stages.children.push_back(leaf(TemporalOps::RhsReady,TemporalOps::Provider));
        for (const auto& p:participants) stages.children.push_back(leaf(p.advance,p.provider,&p));
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
    stage.equationCalls=calls;
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
                                  const std::vector<CompiledMathRef>& calls,
                                  const std::vector<TemporalParticipant>& participants,
                                  const std::vector<SolvePlanNode>& prefix) const override {
        return explicitFragment(body,compiled,calls,participants,prefix);
    }
};

class PartitionedFixedEulerMethod final : public ITemporalMethod {
public:
    FDM::TimeRecipeId id() const override{return FDM::TimeRecipeId::PartitionedFixedEuler;}
    CompiledTimeRecipe compile(FDM::TimeRecipe selected) const override {
        return CompiledTimeRecipe::fromMethod(selected,{{{1.0,0.0,1.0}}},1,ExplicitStageBackend::ForwardEuler);
    }
    SolvePlanNode compileFragment(const CompiledTimeRecipe&,const CompiledEquationCall&,
            const std::vector<CompiledMathRef>&,const std::vector<TemporalParticipant>&,const std::vector<SolvePlanNode>&) const override {
        throw std::runtime_error("partitionedFixedEuler has no RK Stage treatment; declare fixed-time Step calls and HOW loops");
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
                                  const std::vector<CompiledMathRef>& calls,
                                  const std::vector<TemporalParticipant>& participants,
                                  const std::vector<SolvePlanNode>& prefix) const override {
        return explicitFragment(body,compiled,calls,participants,prefix);
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
                                  const std::vector<CompiledMathRef>& calls,
                                  const std::vector<TemporalParticipant>& participants,
                                  const std::vector<SolvePlanNode>& prefix) const override {
        return explicitFragment(body,compiled,calls,participants,prefix);
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
        // The packed viscous kernel reads the model-owned eddy-viscosity
        // array through its transport port. The equation may name that array
        // differently in STATE.
        for(const auto& coefficient:system.state.symbols())
            if(coefficient.storageBinding==StorageBinding::ProviderDistributed
                && coefficient.storageKey=="mu_t" && coefficient.componentOffset==0
                && coefficient.components==1)
                result.coefficientUses.push_back({coefficient.id,{StateVersion::Frozen,StateVersion::Stage}});
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
        call.borrowedViews.push_back({"deltaP",TargetKind::Correction});
    }
    if (method=="ConservativePressurePublication") {
        call.borrowedViews.push_back({"p",TargetKind::Working});
    }
    call.commitRequires={"pressure.published"};
    const bool predictor=method=="ConservativePressureMomentum";
    if (predictor) {
        call.workspaceInvalidates={faceFluxPublication("pressureFaceFlux"),"pressure.predictor","pressureCorrection","preparedPressure","derivedFluxValidity","pressure.published","rAU","faceResponse","correctedFlux"};
        if (conservative) {
            call.stateUses={{"rho"},{"U"},{"p"}};
            call.stateEffects={{"rho"},{"rhoU"},{"rhoE"}};
        }
    } else if (method=="RhieChow") {
        call.workspaceRequires={"rAU"};
        call.stateUses={{"U",StateVersion::Current,false,"",false,TargetKind::Working},{"p"}};
        call.workspaceInvalidates={faceFluxPublication("pressureFaceFlux"),"pressure.published","pressureCorrection","correctedFlux","faceResponse"};
    } else if (method=="PressureCorrection") {
        call.stateUses={{"faceResponse",StateVersion::Frozen,false,"faceResponse"},
            {"correctedFlux",StateVersion::Frozen,false,"correctedFlux"}};
        call.workspaceRequires.push_back("pressure.predictor");
    } else if (method=="ConservativeVelocityCorrection") {
        call.stateEffects={{"rhoU"}};
        call.workspaceInvalidates={"pressure.published"};
    } else if (method=="VelocityCorrection") {
        call.workspaceRequires.push_back("pressure.predictor");
        call.stateUses={{"rAU",StateVersion::Frozen,false,"rAU"}};
        call.workspaceInvalidates={"pressure.published",stateViewPublication("U",StateViewKind::Working)};
        call.stateEffects={{"U",true,true}};
    } else if (method=="PressureUpdate" || method=="ConservativePressureUpdate") {
        call.workspaceInvalidates={"pressure.published"};
        if (!conservative) call.stateEffects={{"p",true,true}};
    } else if (method=="FluxCorrection" || method=="ConservativePressurePublication") {
        if(!conservative)call.stateUses={{"phi",StateVersion::Frozen,false,"correctedFlux"},
            {"faceResponse",StateVersion::Frozen,false,"faceResponse"}};
        call.workspaceInvalidates={"pressureCorrection"};
        if (!conservative) {
            call.workspaceInvalidates.push_back("correctedFlux");
            call.workspaceProvides.push_back("correctedFlux");
        }
        call.workspaceProvides.push_back("pressure.published");
        call.stepWorkspaces.push_back("pressure.published");
        if (!conservative) {
            call.stateEffects={{"phi",true,true}};
            call.faceFlux=CompiledFaceFlux{"phi","U","pressureFaceFlux",TargetKind::Physical,OpIds::PressureCorrectionCommit,{"faceResponse","deltaP"},true};
            call.workspaceProvides.push_back(faceFluxPublication(call.faceFlux->storage));
            call.stepWorkspaces.push_back(faceFluxPublication(call.faceFlux->storage));
        } else call.stateEffects={{"rhoE",true,false}};
    } else if (method=="FixedTimeRelaxation" || method=="ConservativeFixedTimeRelaxation") {
        call.workspaceRequires.push_back("pressure.published");
        call.workspaceInvalidates={"pressure.published",faceFluxPublication("pressureFaceFlux")};
        if (!conservative) call.stateEffects={{"U",true,true},{"p",true,true}};
    } else if (method=="FluxConsistency" || method=="ConservativeFluxConsistency") {
        call.workspaceProvides.push_back("pressure.published");call.stepWorkspaces.push_back("pressure.published");
        if (!conservative) {
            call.stateEffects={{"phi",true,true}};
            call.faceFlux=CompiledFaceFlux{"phi","U","pressureFaceFlux",TargetKind::Physical,OpIds::PressureFluxConsistencyRestore,{"U","p","rAU"},true};
            call.workspaceProvides.push_back(faceFluxPublication(call.faceFlux->storage));
            call.stepWorkspaces.push_back(faceFluxPublication(call.faceFlux->storage));
        }
    } else if (method=="ResidualConvergence" || method=="ConservativeResidualConvergence") {
        call.workspaceRequires.push_back("pressure.published");
    }
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
    if(method=="RhieChow")expected="phi=surfaceFlux(U)";
    else if(method=="PressureCorrection")expected="neg(div(mul(faceResponse,grad(deltaP))))=neg(faceVolumeDivergence(correctedFlux))";
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
                   std::vector<std::string> requirements)
        : id_(method),output_(output),contract_(contract),leaves_(std::move(leaves)),
          needs_(std::move(needs)),provides_(std::move(provides)),
          requirements_(std::move(requirements)) {}
    std::string_view id() const { return id_; }
    std::string_view runtimeProvider() const {
        return (id_=="FluxCorrection" || id_=="RhieChow") ? "flow.rhie-chow" : "flow.pressure-operators";
    }
    std::string_view operationProvider(const OpId& operation) const {
        // Flux publication uses the pressure backend's existing state authority.
        return (operation==OpIds::FluxCorrect || operation==OpIds::FaceFluxPredict) ? "flow.rhie-chow" : "flow.pressure-operators";
    }
    CompiledEquationCall compile(const ExecutableEquationSystem& system,
                                const EquationCall& step,
                                const NumericalBinding& binding) const {
        requireOutput(system,step);
        const bool supportedTarget=step.target.kind==(id_=="RhieChow" ? TargetKind::Working : id_=="PressureCorrection" ? TargetKind::Correction
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
        if (!binding.inputs.empty())
            throw std::runtime_error("EquationMethod '"+std::string(id_)
                +"' has missing or unexpected source-math inputs.");
        CompiledEquationCall result;
        result.source=step;result.backendProvider=runtimeProvider();
        result.target=bindTarget(system,step);
        result.equationMethod=std::string(id_);
        result.target=bindTarget(system,step,provides_.empty() ? "" : provides_.back());
        result.target.viewOwner=StateViewOwner::NumericalProvider;
        result.temporalCapabilities={FDM::TimeRecipeId::ForwardEuler,FDM::TimeRecipeId::PartitionedFixedEuler};
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
        if(id_=="PressureCorrection") {
            result.stepInputs={{OpIds::PressurePrepare,"flow.pressure-operators",false}};
            result.physicalPublication={OpIds::PressureStepCommit,"flow.pressure-operators"};
            result.operations=pressureOperations();
        }
        if(id_=="RhieChow") {
            if(binding.response.empty())throw std::runtime_error("RhieChow needs explicit WHICH response occurrence: "+step.occurrence);
            result.operations.push_back({OpIds::FaceFluxPredict,"construct predicted volume face flux",OperationStage::FluxCorrect,
                {OperationCapability::FluxCorrection},{OriginKind::Generated,"RhieChow"}});
            result.commitRequires.clear(); // A prediction does not require or publish a pressure correction.
            result.workspaceRequires={"cell.response:"+binding.response};
            result.workspaceProvides.push_back("rAU");result.workspaceProvides.push_back("pressure.predictor");
            result.stepWorkspaces=result.workspaceProvides;
            const auto& flux=system.state.at(step.target.symbol);
            Math::requireCompatible(flux.signature,Math::signature(0,{0,3,-1,0,0,0,0}),"RhieChow volume flux");
            if(flux.location!=VariableLocation::EulerianFace || flux.ownership!=OwnershipKind::CanonicalFace || flux.storageKey.empty())
                throw std::runtime_error("RhieChow requires declared canonical scalar face backing: "+step.occurrence);
            if(!step.inputs.empty()) {
                const auto velocity=std::find_if(step.inputs.begin(),step.inputs.end(),[](const auto& use){return use.symbol=="U";});
                if(velocity!=step.inputs.end() && velocity->kind==TargetKind::Physical)result.stateUses.front().kind=TargetKind::Physical;
            }
            if((!step.inputs.empty() && (step.inputs.size()!=result.stateUses.size()
                || !std::equal(step.inputs.begin(),step.inputs.end(),result.stateUses.begin(),[](const auto& a,const auto& b){
                    return a.symbol==b.symbol && a.version==b.version && a.kind==b.kind && a.halo==b.halo && a.snapshot==b.snapshot && a.perIteration==b.perIteration;
                }))))
                throw std::runtime_error("RhieChow requires published Working U/current p and frozen momentum response; no alternate version or inner-stencil override is implemented: "+step.occurrence);
            result.target.workspace="predictedFaceFlux";
            result.target.workspaceLocation=VariableLocation::EulerianFace;
            result.target.workspaceOwnership=OwnershipKind::CanonicalFace;
            result.faceFlux=CompiledFaceFlux{flux.id,"U","predictedFaceFlux",TargetKind::Working,OpIds::FaceFluxPredict,{"U.Working","p","rAU","faceGeometry"},false};
            result.faceFlux->velocityView=result.stateUses.front().kind;
            result.faceFlux->response=binding.response;
            result.faceFlux->numericalInputs.front()=targetText({"U",result.stateUses.front().kind});
            result.spatialTerms=true;
            for(const auto& selected:binding.operators)
                if(selected.provider!="RhieChow")throw std::runtime_error("surfaceFlux occurrence="+selected.occurrence+" has no implementation "+selected.provider+" in selected RhieChow method");
            Math::requireCompatible(system.state.at("U").signature,PressureMathematics::signature("U"),"RhieChow velocity");
            Math::requireCompatible(system.state.at("p").signature,PressureMathematics::signature("p"),"RhieChow pressure");
            result.operatorBindings={"RhieChow"};
        }
        result.fragment.children=leaves_;
        for(auto& leaf:result.fragment.children)leaf.provider=operationProvider(leaf.operation);
        result.fragment.children.back().equationCalls=result.calls;
        if(id_=="PressureUpdate") {
            result.backendOperation=result.fragment.children.back().operation;
            result.residualCapabilities={ResidualKind::StateChange};
        }
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
            result.stepInputs={{OpIds::PressurePrepare,std::string(runtimeProvider()),true}};
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
        result.temporalCapabilities={FDM::TimeRecipeId::ForwardEuler,FDM::TimeRecipeId::PartitionedFixedEuler};
        result.backendOperation=operation_;
        if(id_=="ResidualConvergence") {
            result.residualCapabilities={ResidualKind::DiscreteEquation};result.residualTarget="p";result.residualComponents=2;
        }
        pressureLifecycleContract(result,false);
        if(id_=="FixedTimeRelaxation") {
            result.stepInputs={{OpIds::PressureStepBegin,"flow.pressure-operators"}};
            result.iterationInputs={{OpIds::PressureIterationBegin,"flow.pressure-operators"}};
            result.iterationOutputs={{OpIds::PressureIterationEnd,"flow.pressure-operators"}};
            result.iterationSignalRequired=true;
        }
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

/// @brief 直接方法的选择同时确定本地编译和空间 lowering；无第二套身份匹配。
struct DirectMethod {
    std::string_view name;
    decltype(&compileCellEquation) compile;
    decltype(&lowerCellEquation) lower=nullptr;
};
const DirectMethod* directMethod(std::string_view name) {
    static const DirectMethod methods[]={
        {"LocalAffine",compileLocalAffine},
        {"SSTCentral2",compileSstClosure},
        {"CellExplicit",compileCellEquation,lowerCellEquation},
        {"ConservativeResidual",[](const auto& system,const auto& call,const auto& binding) {
            auto result=compileConservativeResidual(system,call,binding);
            result.backendProvider="flow.conservative";
            return result;
        },[](const auto& system,auto& call,const auto& recipes,const auto& catalog,const auto& bindings,auto& numerical) {
            FusedSpatial::lower(system,call,recipes,catalog,bindings,numerical);
        }}
    };
    const auto selected=std::find_if(std::begin(methods),std::end(methods),[&](const auto& method){return method.name==name;});
    return selected==std::end(methods)?nullptr:selected;
}

} // namespace

TemporalMethodRegistry builtinTemporalMethods() {
    static const ForwardEulerMethod euler;
    static const PartitionedFixedEulerMethod fixed;
    static const SSPRK3Method ssp;
    static const ClassicalRK4Method rk4;
    TemporalMethodRegistry result;
    result.add(euler);result.add(fixed);
    result.add(ssp);
    result.add(rk4);
    return result;
}

std::optional<CompiledEquationCall> compileDirectEquation(const ExecutableEquationSystem& system,
        const EquationCall& call,const NumericalBinding& binding) {
    if(binding.method=="CellExplicitEuler")throw std::runtime_error("Retired method CellExplicitEuler; select CellExplicit and an explicit Step/Stage scope");
    if(const auto* method=directMethod(binding.method))return method->compile(system,call,binding);
    if(auto ordinary=compileScalarEquation(system,call,binding))return ordinary;
    static const PressureMethod prediction("RhieChow","phi",PressureMethod::OutputContract::Assignment,
        {methodLeaf(PlanNodeKind::Update,"RhieChow.predict","predict oriented face volume flux",OpIds::FaceFluxPredict)},
        {"pressure.predictor","rAU"},{"predictedFaceFlux","faceResponse","correctedFlux"},
        {"surfaceFlux(U): oriented integral U dot n dS; SI m^3/s",
         "Working U, current p, frozen rAU, face geometry and original boundary closure",
         "canonical face COPY; requires unchanged ready momentum inputs in the same scope"});
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
    for(const auto* method:{&prediction,&correction,&pressureUpdate,&velocityCorrection,&fluxCorrection})
        if(binding.method==method->id())return method->compile(system,call,binding);
    for(const auto* method:{&relax,&restore,&convergence})
        if(binding.method==method->id())return method->compile(system,call,binding);
    return std::nullopt;
}
bool lowerDirectSpatial(const ExecutableEquationSystem& system,CompiledEquationCall& call,
        const FDM::NumericalRecipeSet& recipes,const TermProviderCatalog& catalog,
        const std::vector<FormulaOperatorBinding>& bindings,CompiledNumericalSystem& result) {
    if(call.equationMethod=="RhieChow") {
        if(result.facePrediction)throw std::runtime_error("Multiple face prediction owners require independently bound face resources");
        const auto& equation=system.registry.at(call.source.equation);
        const auto& surface=equation.rhs;
        const auto occurrence=surface.occurrence().empty()?"rhs":surface.occurrence();
        const auto* selected=selectFormulaBinding(equation,surface,occurrence,bindings);
        if(selected && selected->provider!="RhieChow")
            throw std::runtime_error("equation="+equation.id+" occurrence="+occurrence+" target="+call.target.symbol
                +": selected surfaceFlux implementation="+selected->provider+" is not available in RhieChow");
        result.operators.emplace_back(equation.id,call.source.equation,call.target.symbol,call.equationMethod,
            occurrence,Math::Op::SurfaceFlux,call.faceFlux->velocity,"p",CompiledSpatialBinding::Side::Right,
            std::nullopt,"RhieChow","flow.rhie-chow",true);
        result.facePrediction=call.faceFlux;
        return true;
    }
    const auto* method=directMethod(call.equationMethod);
    if(!method || !method->lower)return false;
    method->lower(system,call,recipes,catalog,bindings,result);
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
    validateScalarExecution(equations,numerics,plan);
    if(!std::isfinite(numerics.dt.maxDeltaT) || numerics.dt.maxDeltaT==std::numeric_limits<double>::max()) {
        for(const auto& term:numerics.operators)
            if(!term.coefficientState.empty() && term.primitiveSpatial==&Discretization::CellEquation::stressArithmeticCentral2)
                throw std::runtime_error("Variable stress needs explicit maxDeltaT before coefficient publication: "+term.callOccurrence+" / "+term.occurrence);
    }
    if(numerics.time.recipe.id()==FDM::TimeRecipeId::PartitionedFixedEuler)
        for(const auto& call:plan.compiledProgram.steps)if(call.temporalResidual || call.stageLocal || call.legacyProvider)
            throw std::runtime_error("partitionedFixedEuler has no Stage/legacy method treatment: "+call.source.occurrence);
    const bool cells=std::any_of(plan.compiledProgram.steps.begin(),plan.compiledProgram.steps.end(),[](const auto& call){return cellEquation(call)!=nullptr && !call.temporalResidual;});
    if(cells && !plan.compiledProgram.temporalParticipants.empty())throw std::runtime_error("Cell Step and recipe-stage integration require a common participant treatment; not implemented");
    if(cells && numerics.time.recipe.id()!=FDM::TimeRecipeId::ForwardEuler && numerics.time.recipe.id()!=FDM::TimeRecipeId::PartitionedFixedEuler)
        throw std::runtime_error("Cell Step requires forwardEuler or partitionedFixedEuler; SSPRK3/classicalRK4 require explicit Stage participation, not a Step response");
    if(!cells && plan.compiledProgram.temporalParticipants.empty()
        && std::any_of(numerics.operators.begin(),numerics.operators.end(),[](const auto& term){return bool(term.fieldSource.bind);}))
        throw std::runtime_error("Field source has no selected numerical reader binding");
    if(numerics.facePrediction) {
        const auto& face=*numerics.facePrediction;
        const auto at=std::find_if(plan.compiledProgram.steps.begin(),plan.compiledProgram.steps.end(),[&](const auto& call){return call.source.occurrence==face.response;});
        const auto* response=at==plan.compiledProgram.steps.end()?nullptr:cellEquation(*at);
        if(!response || at->temporalResidual || at->target.symbol!=face.velocity || response->potential!=face.numericalInputs.at(1)
            || response->gradients.size()!=1 || !response->gradients.front().kernel || !(response->capacityValue>0))
            throw std::runtime_error("RhieChow response="+face.response+" needs a Step Euler response, not an RK participant, with actual negative scalar-gradient/positive-capacity response");
        if(face.velocityView==TargetKind::Working && at->target.kind!=TargetKind::Working
            && std::none_of(at->viewOutputs.begin(),at->viewOutputs.end(),[](const auto& view){return view.kind==StateViewKind::Working;}))
            throw std::runtime_error("RhieChow needs the declared response Working output: "+face.response);
    }
    if(numerics.facePrediction && !numerics.pressureOperator) {
        const auto noPredictedAccept=[&](const auto& self,const ExecutionScope& scope)->void {
            for(const auto& target:scope.commitTargets)if(target.symbol==numerics.facePrediction->symbol)
                throw std::runtime_error("Predicted face Working output is not a corrected physical flux; an explicit implemented correction is required: "+scope.id);
            for(const auto& child:scope.children)self(self,child);
        };noPredictedAccept(noPredictedAccept,plan.compiledProgram.root);
    }
    if(numerics.pressureOperator) {
        const auto reason=validatePressureProgram(equations,numerics,plan);
        if(!reason.empty())throw std::runtime_error(reason);
    }
    if(cells) {
        const bool fixedPressure=std::any_of(plan.compiledProgram.steps.begin(),plan.compiledProgram.steps.end(),[](const auto& call) {
            return std::any_of(call.stepInputs.begin(),call.stepInputs.end(),[](const auto& input) {
                return input.operation==OpIds::PressureStepBegin;
            });
        });
        std::set<std::string> targets;
        for(const auto& call:plan.compiledProgram.steps)if(cellEquation(call) && !call.temporalResidual) {
            if(call.target.kind==TargetKind::Working && (!numerics.pressureOperator || !numerics.facePrediction || numerics.facePrediction->response!=call.source.occurrence)) {
                const auto accepts=[&](const auto& self,const ExecutionScope& scope)->bool {
                    for(std::size_t n=0;n<scope.commitSources.size();++n)if(scope.commitSources[n].symbol==call.target.symbol && scope.commitSources[n].kind==TargetKind::Working
                        && scope.commitTargets.at(n).symbol==call.target.symbol)return true;
                    return std::any_of(scope.children.begin(),scope.children.end(),[&](const auto& child){return self(self,child);});
                };
                if(!accepts(accepts,plan.compiledProgram.root))throw std::runtime_error("Working cell result needs an explicit physical variable Commit: "+call.source.occurrence);
            }
            if(!targets.insert(equations.state.at(call.target.symbol).storageKey).second)
                throw std::runtime_error("Repeated full-dt cell update to one storage needs an explicit fixed-time implementation: "+call.source.occurrence);
        }
        const auto repeated=[&](const auto& self,const ExecutionScope& scope,bool loop)->void {
            loop=loop || ((scope.kind==ExecutionKind::Loop || scope.kind==ExecutionKind::Sequence) && scope.repetitions!=1);
            const bool fixedResponse=fixedPressure && numerics.facePrediction && numerics.facePrediction->response==scope.step.occurrence;
            if(loop && !fixedResponse && scope.kind==ExecutionKind::EquationCall
                && std::any_of(plan.compiledProgram.steps.begin(),plan.compiledProgram.steps.end(),[&](const auto& call){return call.source.occurrence==scope.step.occurrence && cellEquation(call);}))
                throw std::runtime_error("Cell full-dt update cannot repeat without a fixed-time base contract: "+scope.step.occurrence);
            for(const auto& child:scope.children)self(self,child,loop);
        };repeated(repeated,plan.compiledProgram.root,false);
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
    result.add(linearBlockMethod());
    addEulerianMethods(result);
    addImmersedMethods(result);
    addLevelSetMethods(result);
    addMultiphaseMethods(result);
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
