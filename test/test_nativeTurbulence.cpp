#include "models/turbulence/SF_turbulenceSystemContribution.h"
#include "models/turbulence/SF_turbulence.h"
#include "models/physics/fluidStateModel/SF_factory.h"
#include "solver/system/SF_systemBuilder.h"
#include "solver/system/SF_methodObjects.h"
#include "solver/system/SF_builtinState.h"
#include "solver/system/SF_pressureCoupling.h"
#include "solver/system/SF_stateRealizer.h"
#include "solver/run/SF_planExecutor.h"
#include "core/system/SF_operationIds.h"
#include <algorithm>
#include <stdexcept>
using namespace SF::System;
namespace {
void require(bool value,const char* message) { if (!value) throw std::runtime_error(message); }
ExecutableEquationSystem equations(const SystemContribution& contribution) {
    ExecutableEquationSystem system;
    for (auto variable:{BuiltinState::Density,BuiltinState::Momentum,BuiltinState::TotalEnergy,BuiltinState::Velocity}) system.state.add(builtinState(variable));
    for (const auto& state:contribution.states) system.state.add(state);
    for (const auto& equation:contribution.registeredEquations) system.registry.add(equation);
    return system;
}
}
int main() {
    for (const auto* model:{"kEpsilon","kOmegaSST"}) {
        const std::string second=std::string(model)=="kEpsilon" ? "epsilon" : "omega";
        SystemContribution contribution;SF::Turbulence::contribute(contribution,{model,{},false});
        require(contribution.legacyExecution.empty() && contribution.legacyEquations.empty(),"Single-fluid RAS retained legacy authority.");
        require(contribution.registeredEquations.size()==2 && contribution.execution.size()==2
            && contribution.numerics.size()==2 && contribution.states.size()==3,"Missing native RAS channels.");
        auto system=equations(contribution);
        require(system.state.at("mu_t").role==StateRole::Derived && !system.state.isSolution("k"),
            "Model state was confused with user primary selection.");
        for (const auto& equation:system.registry.entries()) {
            require(equation.lhs.name=="ddt","Frozen RAS transient math changed.");
            const auto hasDiv=[](const auto& self,const FormulaExpr& expr)->bool {
                if (expr.kind==FormulaExpr::Kind::Operator && expr.name=="div") return true;
                for (const auto& child:expr.arguments) if (self(self,child)) return true;
                return false;
            };
            require(!hasDiv(hasDiv,equation.rhs),"WHAT falsely claims convection in frozen RAS kernel.");
        }
        ExecutionProgram program;program.root.children=contribution.execution;
        const auto compile=[&](const ExecutionProgram& p,const ExecutableEquationSystem& e) {
            return compileExecutionProgram(e,p,contribution.numerics,builtinProviders());
        };
        const auto compiled=compile(program,system);
        require(compiled.steps.size()==2 && compiled.loweredRoot.children.size()==1
            && compiled.loweredRoot.children[0].operation==OpIds::TurbulenceAdvance
            && compiled.loweredRoot.children[0].provider=="flow.turbulence"
            && compiled.loweredRoot.children[0].equationCalls.size()==2,"RAS pair did not fuse exactly once with provenance.");
        for (const auto& view:compiled.stateViews)
            require(view.owner==StateViewOwner::NumericalProvider,"RAS STATE falsely claims PhysicalState ownership.");
        const auto fails=[&](ExecutionProgram p,ExecutableEquationSystem e,std::vector<NumericalBinding> bindings) {
            bool rejected=false;
            try { (void)compileExecutionProgram(e,p,bindings,builtinProviders()); }
            catch (const std::runtime_error&) { rejected=true; }
            require(rejected,"Invalid RAS fusion contract was accepted.");
        };
        auto missing=program;missing.root.children.pop_back();
        fails(missing,system,{contribution.numerics[0]});
        auto duplicate=program;duplicate.root.children.push_back(program.root.children[0]);
        fails(duplicate,system,contribution.numerics);
        auto wrong=program;wrong.root.children[0].step.target.symbol="mu_t";
        fails(wrong,system,contribution.numerics);
        auto separated=program;ExecutionScope barrier;barrier.kind=ExecutionKind::Sequence;
        separated.root.children[1].order=3;barrier.order=2;separated.root.children.insert(separated.root.children.begin()+1,barrier);
        fails(separated,system,contribution.numerics);
        auto incompatible=system;StateSymbol extra=system.state.at(second);extra.id=second=="omega" ? "epsilon" : "omega";incompatible.state.add(extra);
        fails(program,incompatible,contribution.numerics);
        auto changed=system;auto formula=system.registry.at("k");formula.rhs=FormulaExpr::constantValue(0);changed.registry.replace(formula);
        fails(program,changed,contribution.numerics);

        SF::FDM::SolverConfig config;config.numerics.timeRecipe=SF::FDM::builtInTimeRecipe(SF::FDM::TimeRecipeId::ClassicalRK4);
        config.numerics.recipes.time=config.numerics.timeRecipe;
        BuildRequest request;request.singleFluidPreset=SingleFluidPresetSpec{false};
        request.composition.stateDeclared=true;request.composition.solutionVariables={"rho","rhoU","rhoE"};
        request.modelContributions={contribution};
        const auto resolved=build(config,request);
        require(resolved.solvePlan.compiledProgram.steps.size()==5 && resolved.solvePlan.compiledProgram.hasTemporalRoot,
            "Native RAS and flow did not compose.");
        SF::Run::OpRegistry operations;std::vector<std::string> order;int advances=0,stages=0;
        const auto bind=[&](const auto& self,const SolvePlanNode& node)->void {
            if (!node.operation.empty()) operations.bind(node.operation,node.provider,[&,op=node.operation] {
                order.push_back(op);advances+=op==OpIds::TurbulenceAdvance;stages+=op==OpIds::ExplicitStageExecute;
            });
            for (const auto& child:node.children) self(self,child);
        };
        bind(bind,resolved.solvePlan.root);SF::Run::PlanExecutor::execute(resolved.solvePlan,operations);
        require(advances==1 && stages==4 && order[2]==OpIds::TurbulenceAdvance && order[3]==OpIds::FlowStepBegin,
            "RAS was double-advanced or moved relative to flow snapshot.");
        auto distributed=request;distributed.parallel=true;
        const auto unsupported=build(config,distributed);
        const auto scope=std::find_if(unsupported.runtime.requirements.begin(),unsupported.runtime.requirements.end(),
            [](const auto& requirement) { return requirement.name=="SingleFluidTurbulenceScope"; });
        require(scope!=unsupported.runtime.requirements.end() && !scope->available,"RAS MPI capability was silently expanded.");
        auto pressure=request;pressure.singleFluidPreset.reset();
        pressure.composition.declared=true;pressure.composition.solutionVariables={"U","p"};
        pressure.composition.equations={"Continuity","Momentum"};pressure.composition.algorithm="PISO";
        pressure.composition.thermoDynamics.equationOfState="rhoConst";
        pressure.composition.thermoDynamics.transport="const";
        pressure.composition.thermoDynamics.constantDensity=1;
        config.numerics.timeRecipe=SF::FDM::builtInTimeRecipe(SF::FDM::TimeRecipeId::ForwardEuler);
        config.numerics.recipes.time=config.numerics.timeRecipe;
        pressure.coupling=couplingRequestFrom(config.pressure.coupling,true);
        pressure.coupling->presetKind=SF::FDM::PressureCouplingPreset::PISO;
        const auto pressureSystem=build(config,pressure);
        const auto rejected=std::find_if(pressureSystem.runtime.operationBindings.begin(),pressureSystem.runtime.operationBindings.end(),
            [](const auto& binding) { return binding.operation==OpIds::TurbulenceAdvance; });
        require(rejected!=pressureSystem.runtime.operationBindings.end() && rejected->status==BindingStatus::Unsupported,
            "Pressure plus RAS was falsely advertised as executable.");
    }
    for (const auto* model:{"DNS","Smagorinsky"}) {
        SystemContribution contribution;SF::Turbulence::contribute(contribution,{model,{},false});
        require(std::none_of(contribution.states.begin(),contribution.states.end(),[](const auto& state) {
            return state.role==StateRole::Transported;
        }),"Closure-only model acquired transported state.");
        require(std::none_of(contribution.numerics.begin(),contribution.numerics.end(),[](const auto& method) {
            return method.method=="TurbulenceTransport";
        }),"Closure-only model acquired RAS transport method.");
    }
    SystemContribution eulerian;SF::Turbulence::contribute(eulerian,{"kOmegaSST",{"water"},true,{"water","air"}});
    require(eulerian.legacyExecution.empty() && eulerian.legacyEquations.empty() && eulerian.execution.size()==3,"Eulerian native model channels are incomplete.");

    // Real provider-owned backing: STATE views alias ScalarFields, including writes.
    SF::Field field;field.setup(3,3,1,1,5);
    field.setStateModel(SF::Physics::FluidStateModel::makeSingleFluidPerfectGas(1.4,287.05,0.0,0.72));
    for (int cell=0;cell<field.TotalSize();++cell) {
        int i,j,k;field.getIJK(cell,i,j,k);field(i,j,k,0)=1;field(i,j,k,4)=2.5;
        field.Jac(i,j,k)=1;field.XiX(i,j,k)=1;field.EtY(i,j,k)=1;field.ZeZ(i,j,k)=1;
    }
    SF::FDM::TurbulenceConfig config;config.enabled=true;config.family=SF::FDM::TurbulenceFamily::RAS;
    config.model=SF::FDM::TurbulenceModelKind::kEpsilon;
    SF::Turbulence::Manager manager(config);require(manager.initialize(field),"Manager initialization failed.");
    SF::State::StateBundle bundle;bundle.patches={&field};manager.registerDistributed(bundle.distributed,0,field);
    SystemContribution contribution;SF::Turbulence::contribute(contribution,{"kEpsilon",{},false});
    const auto system=equations(contribution);ExecutionProgram program;program.root.children=contribution.execution;
    const auto compiled=compileExecutionProgram(system,program,contribution.numerics,builtinProviders());
    const auto realized=realizeState(system.state,RuntimeRequirements{},bundle,compiled.stateViews);
    const auto& view=realized.view("k",StateViewKind::Physical);
    view.fields[0]->write(0,0,0.25);
    require(manager.scalarFields().values(SF::Turbulence::ScalarSlot::K)[0]==0.25
        && realized.at("k").fields[0]==view.fields[0],"STATE allocated duplicate turbulence backing.");
    std::fill(manager.scalarFields().values(SF::Turbulence::ScalarSlot::K).begin(),
        manager.scalarFields().values(SF::Turbulence::ScalarSlot::K).end(),0.01);
    std::fill(manager.scalarFields().values(SF::Turbulence::ScalarSlot::Epsilon).begin(),
        manager.scalarFields().values(SF::Turbulence::ScalarSlot::Epsilon).end(),0.001);
    std::vector<double> before;
    for (int cell=0;cell<field.TotalSize();++cell) {
        int i,j,k;field.getIJK(cell,i,j,k);
        for (int c=0;c<5;++c) before.push_back(field(i,j,k,c));
    }
    SF::Run::OpRegistry operations;CompiledSolvePlan plan;plan.root=compiled.loweredRoot;
    operations.bind(OpIds::TurbulenceAdvance,"flow.turbulence",[&] { manager.correct(field,0.01); });
    SF::Run::PlanExecutor::execute(plan,operations);
    require(manager.scalarFields().K(1,1,1)!=0.01 && manager.scalarFields().Epsilon(1,1,1)!=0.001
        && manager.scalarFields().EddyMu(1,1,1)>0,"Fused RAS operation did not publish transported/closure arrays.");
    std::size_t n=0;
    for (int cell=0;cell<field.TotalSize();++cell) {
        int i,j,k;field.getIJK(cell,i,j,k);
        for (int c=0;c<5;++c) require(before[n++]==field(i,j,k,c),"RAS directly overwrote flow primary state.");
    }
}
