#include "solver/system/SF_systemBuilder.h"
#include "solver/system/SF_methodObjects.h"
#include "solver/system/SF_pressureCoupling.h"
#include "solver/system/SF_stateRealizer.h"
#include "solver/run/SF_planExecutor.h"
#include "SF_phaseState.h"
#include "models/turbulence/SF_turbulenceSystemContribution.h"
#include "core/system/SF_operationIds.h"
#include <algorithm>
#include <map>
#include <stdexcept>
using namespace SF::System;
namespace {
void require(bool value,const char* message) { if (!value) throw std::runtime_error(message); }
void calls(ExecutionScope& scope,const std::string& id,const std::function<void(ExecutionScope&)>& edit) {
    if (scope.kind==ExecutionKind::EquationCall && scope.step.equation==id) edit(scope);
    for (auto& child:scope.children) calls(child,id,edit);
}
}
int main() {
    SF::FDM::SolverConfig config;config.numerics.maxDeltaT=0.001;
    config.pressure.coupling.preset=SF::FDM::PressureCouplingPreset::PIMPLE;
    config.pressure.coupling.outerCorrectors=2;config.pressure.coupling.pressureCorrectors=3;
    config.pressure.coupling.nonOrthogonalCorrectors=1;
    BuildRequest request;request.templateOrigin=PhysicsTemplateKind::EulerianEulerian;
    request.phaseNames={"water","air","other"};request.referencePhase="water";
    request.coupling=couplingRequestFrom(config.pressure.coupling,true);
    const auto resolved=build(config,request);
    require(resolved.executionPolicies.empty() && resolved.solvePlan.sourceProgram.legacyEntries.empty(),"Eulerian retained a second schedule authority.");
    require(resolved.executableSystem.state.solutionVariables().empty(),"Phase fields became single-fluid primary selection.");
    require(resolved.solvePlan.compiledProgram.steps.size()==17,"Missing separate phase/correction occurrences.");
    require(!resolved.executableSystem.state.contains("pPrime") && !resolved.executableSystem.state.contains("pressureCorrection"),"Pressure correction became another base state.");
    for (const auto& binding:resolved.runtime.operationBindings)
        require(binding.provider=="flow.eulerian-pressure" && binding.status==BindingStatus::Resolved,"Provider ownership was not frozen.");
    const auto& root=resolved.solvePlan.sourceProgram.root;
    require(root.children.size()==2 && root.children[0].repetitions==2
        && root.children[0].children[2].repetitions==3
        && root.children[0].children[2].children[0].repetitions==2,"Native HOW lost nested loop counts.");
    SF::Run::OpRegistry operations;std::map<std::string,int> counts;std::vector<std::string> order;
    for (const auto& id:resolved.runtime.report.requiredOperations)
        operations.bind(id,"flow.eulerian-pressure",[&,id] { ++counts[id];order.push_back(id); });
    SF::Run::PlanExecutor::execute(resolved.solvePlan,operations);
    require(counts[OpIds::EeDtCompute]==1 && counts[OpIds::EeStepBegin]==1
        && counts[OpIds::EeContinuityAssemble]==2 && counts[OpIds::EeMomentumSolve]==2
        && counts[OpIds::EeEnergySolve]==2 && counts[OpIds::EePressureSolve]==12
        && counts[OpIds::EePhaseCorrect]==12 && counts[OpIds::EeFaceFluxCorrect]==12
        && counts[OpIds::EeStepCommit]==1 && counts[OpIds::EeTimeCommit]==1,"Multi-phase kernels executed once per member instead of group.");
    const auto provenance=[&](const auto& self,const SolvePlanNode& node)->void {
        if (!node.operation.empty()) {
            require(!node.legacyAdapter && node.provider=="flow.eulerian-pressure","Core leaf uses legacy provider discovery.");
            if (node.operation==OpIds::EeContinuityAssemble || node.operation==OpIds::EeMomentumSolve || node.operation==OpIds::EeEnergySolve)
                require(node.equationCalls.size()==3,"Grouped phase provenance disappeared.");
            if (node.operation==OpIds::EePhaseCorrect) require(node.equationCalls.size()==4,"Shared/phase correction provenance disappeared.");
        }
        for (const auto& child:node.children) self(self,child);
    };
    provenance(provenance,resolved.solvePlan.root);
    const auto reject=[&](ExecutionProgram program,ExecutableEquationSystem equations,std::vector<NumericalBinding> bindings) {
        bool rejected=false;try { (void)compileExecutionProgram(equations,program,bindings,builtinProviders()); }
        catch (const std::runtime_error&) { rejected=true; }
        require(rejected,"Invalid Eulerian composition was accepted.");
    };
    auto missing=resolved.solvePlan.sourceProgram;
    missing.root.children[0].children[0].children.pop_back();
    reject(missing,resolved.executableSystem,resolved.numericalSelection.bindings);
    auto duplicate=resolved.solvePlan.sourceProgram;
    duplicate.root.children[0].children[1].children.push_back(duplicate.root.children[0].children[1].children.front());
    reject(duplicate,resolved.executableSystem,resolved.numericalSelection.bindings);
    auto wrong=resolved.solvePlan.sourceProgram;
    calls(wrong.root,"momentum.air",[](auto& node) {node.step.target.symbol="U.air";});
    reject(wrong,resolved.executableSystem,resolved.numericalSelection.bindings);
    auto cross=resolved.solvePlan.sourceProgram;
    auto moved=cross.root.children[0].children[1].children.back();cross.root.children[0].children[1].children.pop_back();
    cross.root.children[0].children.insert(cross.root.children[0].children.begin()+2,moved);
    reject(cross,resolved.executableSystem,resolved.numericalSelection.bindings);
    auto changed=resolved.executableSystem;auto formula=changed.registry.at("momentum.air");formula.rhs=FormulaExpr::constantValue(0);changed.registry.replace(formula);
    reject(resolved.solvePlan.sourceProgram,changed,resolved.numericalSelection.bindings);
    auto badBacking=resolved.executableSystem;badBacking.state=StateRegistry{};
    for (auto symbol:resolved.executableSystem.state.symbols()) {
        if (symbol.id=="momentum.air") symbol.storageKey="phase0.momentum";
        badBacking.state.add(symbol);
    }
    reject(resolved.solvePlan.sourceProgram,badBacking,resolved.numericalSelection.bindings);
    auto noCommit=resolved.solvePlan.sourceProgram;noCommit.root.children.pop_back();
    reject(noCommit,resolved.executableSystem,resolved.numericalSelection.bindings);
    auto reordered=resolved.solvePlan.sourceProgram;
    std::swap(reordered.root.children[0].children[1],reordered.root.children[0].children[3]);
    reject(reordered,resolved.executableSystem,resolved.numericalSelection.bindings);
    auto correction=resolved.solvePlan.sourceProgram;
    calls(correction.root,"E_SHARED_PRESSURE",[](auto& node) {node.step.target.kind=TargetKind::Physical;});
    reject(correction,resolved.executableSystem,resolved.numericalSelection.bindings);
    auto unsupported=request;unsupported.coupling->presetKind=SF::FDM::PressureCouplingPreset::PISO;
    unsupported.coupling->preset="PISO";bool rejected=false;
    try { (void)build(config,unsupported); } catch (const std::runtime_error&) { rejected=true; }
    require(rejected,"Eulerian silently selected PIMPLE for requested PISO.");
    auto legacy=request;SystemContribution turbulence;
    SF::Turbulence::contribute(turbulence,{"kOmegaSST",request.phaseNames,true});
    legacy.modelContributions.push_back(turbulence);
    const auto compatibility=build(config,legacy);
    require(!compatibility.solvePlan.sourceProgram.legacyEntries.empty() && compatibility.executionPolicies.empty()
        && compatibility.runtime.report.status==RuntimeStatus::Unsupported
        && compatibility.runtime.report.reason.find("has no compiled HOW execution step")!=std::string::npos,
        "Unmigrated Eulerian turbulence was silently enabled/disabled or acquired another scheduler.");
    // The real PhaseState arrays and correction scalar remain the sole backing authorities.
    SF::Field field;field.setup(2,2,1,1,5);SF::State::StateBundle bundle;bundle.patches={&field};
    std::vector<SF::Physics::PhaseSystems::PhaseState> phases(3);
    SF::ScalarField pressure,delta;pressure.setupLike(field,"p",101325);delta.setupLike(field,"delta",0);
    std::vector<std::vector<double>> flux(3,std::vector<double>(3*field.TotalSize()));
    bundle.distributed.add(SF::State::scalarView("pressure",0,field,pressure));
    bundle.distributed.add(SF::State::scalarView("pressureCorrection",0,field,delta));
    for (std::size_t i=0;i<phases.size();++i) {
        phases[i].setupLike(field,request.phaseNames[i]);auto& p=phases[i];const auto key="phase"+std::to_string(i)+".";
        for (const auto& item:std::vector<std::pair<std::string,SF::ScalarField*>>{
                {"mass",&p.primary.phaseMass},{"enthalpy",&p.primary.phaseEnthalpy},{"alpha",&p.primitive.alpha},
                {"density",&p.primitive.density},{"primitiveEnthalpy",&p.primitive.enthalpy},{"temperature",&p.primitive.temperature}})
            bundle.distributed.add(SF::State::scalarView(key+item.first,0,field,*item.second));
        for (const auto& item:std::vector<std::pair<std::string,SF::Physics::PhaseSystems::PhaseVectorField*>>{
                {"momentum",&p.primary.momentum},{"velocity",&p.primitive.velocity}})
            bundle.distributed.add(SF::State::scalarComponentsView(key+item.first,0,field,{&(*item.second)[0],&(*item.second)[1],&(*item.second)[2]},0,SF::State::HaloSyncStage::None));
        bundle.distributed.add(SF::State::workspaceView(key+"volumeFaceFlux",0,field,flux[i],3));
    }
    auto realized=realizeState(resolved.executableSystem.state,{},bundle,resolved.solvePlan.compiledProgram.stateViews);
    auto* mass=realized.view("phaseMass.air",StateViewKind::Physical).fields.front();mass->write(0,0,0.42);
    require(phases[1].primary.phaseMass.values()[0]==0.42,"Phase STATE allocated a second authority.");
    auto deltas=bundle.distributed.select("pressureCorrection",SF::State::HaloSyncStage::None);
    realized.bindView("p",StateViewKind::Correction,*deltas.front());
    realized.view("p",StateViewKind::Correction).fields.front()->write(0,0,1.5);
    require(delta.values()[0]==1.5 && pressure.values()[0]==101325,"Correction(p) aliases or overwrites physical pressure incorrectly.");
}
