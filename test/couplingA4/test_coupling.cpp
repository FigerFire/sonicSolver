#include "solver/system/SF_systemBuilder.h"
#include "solver/system/SF_systemValidator.h"
#include "solver/system/SF_systemPrinter.h"
#include "solver/system/SF_methodObjects.h"
#include "solver/system/SF_executionContract.h"
#include "solver/system/SF_formulaCompiler.h"
#include "solver/system/SF_solvePlan.h"
#include "solver/run/SF_planExecutor.h"
#include "core/system/SF_operationIds.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <filesystem>
using namespace SF;
using namespace SF::System;
namespace {
int checks=0;
void require(bool ok,const std::string& why) {++checks;if(!ok)throw std::runtime_error(why);}
template<class F> void rejects(const std::string& text,F action) {
    try {action();} catch(const std::runtime_error& e) {
        require(std::string(e.what()).find(text)!=std::string::npos,e.what());
        std::cout<<"Rejected ["<<text<<"] "<<e.what()<<'\n';return;
    }
    throw std::runtime_error("Expected rejection: "+text);
}
BuildRequest request(bool conservative,FDM::PressureCouplingPreset preset) {
    BuildRequest r;r.composition.stateDeclared=true;
    r.composition.solutionVariables=conservative ? std::vector<std::string>{"rho","rhoU","rhoE"} : std::vector<std::string>{"U","p"};
    if(conservative) r.pressureConstraint=PressureConstraintSpec{};
    else {
        r.composition.declared=true;r.composition.algorithm=FDM::toString(preset);r.composition.equations={"Momentum","Continuity"};
        r.composition.thermoDynamics.equationOfState="rhoConst";
        r.composition.thermoDynamics.transport="const";r.composition.thermoDynamics.constantDensity=1;
    }
    FDM::PressureCouplingConfig c;c.preset=preset;
    c.outerCorrectors=preset==FDM::PressureCouplingPreset::PISO ? 1 : 3;
    c.pressureCorrectors=preset==FDM::PressureCouplingPreset::SIMPLE ? 1 : 2;
    r.coupling=couplingRequestFrom(c,true);return r;
}
void renamed(bool conservative,FDM::PressureCouplingPreset preset) {
    FDM::SolverConfig config;config.numerics.timeRecipe=FDM::builtInTimeRecipe(FDM::TimeRecipeId::ForwardEuler);
    if (!conservative) {
        config.numerics.formulation=FDM::EquationFormulation::PrimitiveDifferential;
        config.numerics.convection=FDM::ConvectionScheme::Upwind1;
        config.numerics.reconstruction=FDM::ReconstructionVariable::Primitive;
        config.numerics.flux=FDM::FluxSplitter::UpwindAdvection;
        config.numerics.recipes.convection=FDM::TermRecipe::primitiveUpwindRecipe();
    }
    auto r=request(conservative,preset);
    const auto original=build(config,r);validate(original);
    SystemContribution authored;authored.recordContribution("test.renamed","explicit pressure references");
    for(const auto& state:original.rawSystem.state.symbols()) authored.addState(state);
    for(auto equation:original.rawSystem.registry.entries()) {
        if(equation.id=="momentum") equation.id="user.vectorBalance";
        authored.addEquation(std::move(equation));
    }

    // This is an explicitly authored input occurrence, with no default node ID or preset provenance.
    auto input=pressureCouplingInputs({"user.vectorBalance",{conservative ? "rhoU" : "U"},"user.seed"},
        conservative ? PressurePredictorKind::ConservativeMomentum : PressurePredictorKind::Velocity);
    input.predictorOccurrence="user.predict";
    std::vector<std::pair<std::string,std::string>> renamedRelations;
    for(auto* relation:{&input.correction,&input.velocityUpdate,&input.pressureUpdate,&input.fluxUpdate,
            &input.publication,&input.relaxation,&input.fluxConsistency,&input.convergence}) {
        if(relation->equation.empty()) continue;
        renamedRelations.push_back({relation->equation,"user."+relation->equation});
        relation->equation="user."+relation->equation;
        relation->occurrence="occurrence."+relation->equation;
    }
    for(auto constraint:original.rawSystem.constraints) {
        for(auto& member:constraint.equations) for(const auto& rename:renamedRelations)
            if(member==rename.first) member=rename.second;
        authored.addConstraint(constraint);
    }
    ExecutionScope predictor;predictor.kind=ExecutionKind::EquationCall;
    predictor.step=input.predictor;predictor.id="arbitrary-node-name";authored.addExecution(predictor);
    for(const auto& balance:input.fusedBalances) {
        ExecutionScope call;call.kind=ExecutionKind::EquationCall;call.step=balance;authored.addExecution(call);
    }
    r.includeDefaultFluidPreset=false;r.userContributions={authored};r.pressureInputs=input;
    const auto system=build(config,r);validate(system);
    require(original.runtime.report.status==RuntimeStatus::Runnable && system.runtime.report.status==RuntimeStatus::Runnable,"Rename changed runnable provider availability");
    require(SolvePlanner::requiredOperations(system.solvePlan)==SolvePlanner::requiredOperations(original.solvePlan),"Rename changed required operations");
    for(const auto& compiled:system.solvePlan.compiledProgram.steps) for(const auto& view:compiled.viewReads)
        require(view.owner==StateViewOwner::NumericalProvider && !view.storage.empty(),"Pressure read view was not frozen to the existing resource");
    const auto& call=system.solvePlan.compiledProgram.steps.front();
    require(call.source.equation=="user.vectorBalance" && call.source.occurrence=="user.predict","Compiler lost explicit identity");
    require(canonicalFormula(system.rawSystem.registry.at("user.vectorBalance"))==canonicalFormula(original.rawSystem.registry.at("momentum")),"Rename changed mathematics");
    require(describe(system).find("transforms user.vectorBalance occurrence=user.seed")!=std::string::npos,"Explain mislabels explicit input as a default");
    require(describe(system).find("user.predict")!=std::string::npos,"Explain lost explicit reference");
    std::vector<std::string> trace;
    Run::OpRegistry ops;
    const auto bind=[&](const auto& self,const SolvePlanNode& node)->void {
        if(!node.operation.empty() && !ops.contains(node.operation)) ops.bind(node.operation,node.provider,[&,id=node.operation](const Run::ExecutionContext& c){
            trace.push_back(id);
            if(id==OpIds::PressureConvergenceEvaluate) c.signals->publish(kPressureOuterConvergedSignal,false);
        });
        for(const auto& child:node.children) self(self,child);
    };
    bind(bind,system.solvePlan.root);Run::PlanExecutor::execute(system.solvePlan,ops);
    const auto renamedTrace=trace;trace.clear();Run::PlanExecutor::execute(original.solvePlan,ops);
    require(trace==renamedTrace,"Renamed source changed actual loop execution");
    auto invalid=r;invalid.pressureInputs->correction.target.kind=TargetKind::Working;
    rejects("semantic target",[&]{build(config,invalid);});
    invalid=r;invalid.pressureInputs->predictor.occurrence="missing";
    rejects("referenced predictor",[&]{build(config,invalid);});
    auto compiled=system.solvePlan.compiledProgram;
    auto& calls=compiled.steps;
    const auto publisher=std::find_if(calls.begin(),calls.end(),[&](const auto& c){return c.equationMethod==(conservative ? "ConservativePressurePublication" : "FluxCorrection");});
    require(publisher!=calls.end(),"Missing real publication");
    publisher->workspaceProvides.erase(std::remove(publisher->workspaceProvides.begin(),publisher->workspaceProvides.end(),"pressure.published"),publisher->workspaceProvides.end());
    rejects(preset==FDM::PressureCouplingPreset::PISO ? "Commit requires unpublished" : "requires unavailable workspace",[&]{validateDataFlow(system.executableSystem.state,compiled);});
    compiled=system.solvePlan.compiledProgram;
    const auto correction=std::find_if(compiled.steps.begin(),compiled.steps.end(),[&](const auto& c){return c.equationMethod==(conservative ? "ConservativePressureCorrection" : "PressureCorrection");});
    correction->workspaceProvides.clear();correction->stepWorkspaces.clear();
    rejects("requires unavailable workspace",[&]{validateDataFlow(system.executableSystem.state,compiled);});
    if(!conservative) {
        compiled=system.solvePlan.compiledProgram;
        compiled.steps[1].stateUses={{"U",StateVersion::Stage}};
        rejects("outside declared StageLoop",[&]{validateDataFlow(system.executableSystem.state,compiled);});
    }
    if(conservative) {
        compiled=system.solvePlan.compiledProgram;
        for(auto& effect:compiled.steps.front().stateEffects) effect.invalidatesLazy=false;
        CompiledEquationCall read;read.source={"observer",{"rhoU"},"observer.current"};read.stateUses={{"U"}};
        ExecutionScope node;node.kind=ExecutionKind::EquationCall;node.step=read.source;
        compiled.root.children.insert(compiled.root.children.begin()+1,node);compiled.steps.push_back(read);
        rejects("Lazy STATE invalidation missing",[&]{validateDataFlow(system.executableSystem.state,compiled);});
    }
    std::ofstream("result/"+std::string(conservative ? "conservative-" : "primitive-")+FDM::toString(preset)+"-explain.txt")<<describe(system);
}
void genericIteration() {
    ExecutableEquationSystem system;
    for(const auto* id:{"a","b"}) {
        StateSymbol s;s.id=s.name=s.storageKey=id;s.storageBinding=StorageBinding::NamedDistributed;
        s.signature=Math::signature(0,{});system.state.add(s);
    }
    const auto a=system.state.at("a").expression(),b=system.state.at("b").expression();
    system.registry.add({"update-a",a,FormulaExpr::constantValue(.5)*(b+FormulaExpr::constantValue(1))});
    system.registry.add({"update-b",b,FormulaExpr::constantValue(.5)*(a+FormulaExpr::constantValue(1))});
    ExecutionProgram how;how.explicitOrder=true;
    ExecutionScope loop;loop.kind=ExecutionKind::Loop;loop.id="arbitrary-iteration";loop.repetitions=12;
    for(auto call:std::vector<EquationCall>{{"update-a",{"a"},"first"},{"update-b",{"b"},"second"}}) {
        ExecutionScope node;node.kind=ExecutionKind::EquationCall;node.step=call;loop.children.push_back(node);
    }
    how.root.children={loop};
    auto compiled=compileExecutionProgram(system,how,{{"update-a","DirectEvaluation"},{"update-b","DirectEvaluation"}},builtinProviders());
    // Bind existing generic value kernels to two scalar resources. No CFD/time loop or pressure callback.
    double av=0,bv=0;FormulaValues values;values.boundReads={ [&](int,int){return av;},[&](int,int){return bv;} };
    FormulaOperatorCatalog catalog;Run::OpRegistry ops;
    for(auto& call:compiled.steps) {
        const auto& equation=system.registry.at(call.source.equation);std::vector<std::string> used;
        auto kernel=compileFormulaValue(equation,equation.rhs,catalog,used,{"a","b"},&system.state);
        const bool first=call.target.symbol=="a";
        ops.bind(call.backendOperation,[&,kernel,first]{(first ? av : bv)=kernel(values,0,0);});
        call.stateUses={{first ? "b" : "a"}};call.stateEffects={{first ? "a" : "b"}};
    }
    validateDataFlow(system.state,compiled);
    CompiledSolvePlan plan;plan.root=compileMethodProgram(compiled);
    Run::PlanExecutor::execute(plan,ops);
    double expectedA=0,expectedB=0;
    for(int i=0;i<12;++i){expectedA=.5*(expectedB+1);expectedB=.5*(expectedA+1);}
    require(av==expectedA && bv==expectedB,"Generic two-equation Loop changed sequential iteration");
    require(std::abs(av-1)<1e-6 && std::abs(bv-1)<1e-6,"Iteration did not approach its algebraic fixed point");
    std::ofstream("result/iteration.txt")<<std::setprecision(17)<<av<<' '<<bv<<'\n';
    compiled.steps[0].workspaceProvides={"iteration.local"};compiled.root.children.push_back(compiled.root.children[0].children.back());
    compiled.root.children.back().step.occurrence="outside";
    auto outside=compiled.steps.back();outside.source=compiled.root.children.back().step;outside.workspaceRequires={"iteration.local"};compiled.steps.push_back(outside);
    rejects("requires unavailable workspace",[&]{validateDataFlow(system.state,compiled);});
    // Invalidating a previously valid lease inside a loop cannot resurrect it on exit.
    CompiledExecutionProgram lease;
    for(const auto* id:{"publish","invalidate","read"}) {
        CompiledEquationCall c;c.source={"update-a",{"a"},id};lease.steps.push_back(c);
        ExecutionScope n;n.kind=ExecutionKind::EquationCall;n.step=c.source;lease.root.children.push_back(n);
    }
    lease.steps[0].workspaceProvides={"lease"};lease.steps[1].workspaceInvalidates={"lease"};lease.steps[2].workspaceRequires={"lease"};
    ExecutionScope scope;scope.kind=ExecutionKind::Loop;scope.id="lease-scope";scope.children={lease.root.children[1]};lease.root.children[1]=scope;
    rejects("requires unavailable workspace",[&]{validateDataFlow(system.state,lease);});
}
}
int main() {
    std::filesystem::create_directories("result");
    for(auto preset:{FDM::PressureCouplingPreset::PISO,FDM::PressureCouplingPreset::SIMPLE,FDM::PressureCouplingPreset::PIMPLE}) renamed(false,preset);
    renamed(true,FDM::PressureCouplingPreset::PISO);
    genericIteration();
    std::cout<<checks<<" A4 checks passed\n";
}
