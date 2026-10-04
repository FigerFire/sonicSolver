#include "models/turbulence/SF_turbulenceSystemContribution.h"
#include "solver/system/SF_systemBuilder.h"
#include "solver/system/SF_eulerianTurbulence.h"
#include "solver/system/SF_eulerianAssembly.h"
#include "solver/system/SF_pressureCoupling.h"
#include "solver/system/SF_stateRealizer.h"
#include "solver/run/SF_planExecutor.h"
#include "core/system/SF_operationIds.h"
#include <algorithm>
#include <map>
#include <stdexcept>
using namespace SF::System;
namespace {
void require(bool value,const char* message) {if (!value) throw std::runtime_error(message);}
void leaves(const SolvePlanNode& node,const std::function<void(const SolvePlanNode&)>& visit) {
    if (!node.operation.empty()) visit(node);
    for (const auto& child:node.children) leaves(child,visit);
}
}
int main() {
    for (const auto* model:{"kOmegaSST","kEpsilon","Smagorinsky"})
        for (const auto& selected:std::vector<std::vector<std::string>>{{"water","air"},{"air"},{"air","water"}}) {
            const bool ras=std::string(model)!="Smagorinsky";
            const auto second=std::string(model)=="kEpsilon" ? "epsilon" : "omega";
            SF::FDM::SolverConfig config;config.numerics.maxDeltaT=.001;
            config.pressure.coupling.preset=SF::FDM::PressureCouplingPreset::PIMPLE;
            config.pressure.coupling.outerCorrectors=2;config.pressure.coupling.pressureCorrectors=2;
            config.pressure.coupling.nonOrthogonalCorrectors=1;
            BuildRequest request;request.templateOrigin=PhysicsTemplateKind::EulerianEulerian;
            request.phaseNames={"water","air"};request.referencePhase="water";
            request.coupling=couplingRequestFrom(config.pressure.coupling,true);
            SystemContribution turbulence;SF::Turbulence::contribute(turbulence,{model,selected,true,request.phaseNames});
            require(turbulence.legacyExecution.empty() && turbulence.legacyEquations.empty(),"Native model retained compatibility entries.");
            request.modelContributions={turbulence};const auto resolved=build(config,request);
            require(resolved.runtime.report.status==RuntimeStatus::Runnable && resolved.solvePlan.sourceProgram.legacyEntries.empty(),"Native Eulerian model is not executable.");
            require(resolved.executableSystem.legacyEquations.empty() && resolved.executableSystem.legacyDefinitions.equations().empty(),"Duplicate legacy math authority.");
            SF::FDM::TurbulenceConfig tc;tc.enabled=true;tc.phaseNames=selected;
            tc.family=ras ? SF::FDM::TurbulenceFamily::RAS : SF::FDM::TurbulenceFamily::LES;
            tc.model=std::string(model)=="kEpsilon" ? SF::FDM::TurbulenceModelKind::kEpsilon
                : ras ? SF::FDM::TurbulenceModelKind::kOmegaSST : SF::FDM::TurbulenceModelKind::Smagorinsky;
            validateEulerianTurbulenceBindings(resolved.solvePlan.compiledProgram,request.phaseNames,tc);
            validateEulerianAssemblyBindings(resolved.solvePlan.compiledProgram,request.phaseNames,0);
            std::map<std::string,int> counts;std::vector<std::string> order;SF::Run::OpRegistry operations;
            leaves(resolved.solvePlan.root,[&](const auto& node) {
                require(!node.legacyAdapter && !node.operation.empty(),"Legacy/unbound native leaf.");
                const bool tur=node.operation==OpIds::EeTurbulencePrepare || node.operation==OpIds::EeTurbulenceSolve;
                require(node.provider==(tur ? "flow.eulerian-turbulence" : "flow.eulerian-pressure"),"Wrong operation ownership.");
                if (tur) require(node.equationCalls.size()==selected.size()*(node.operation==OpIds::EeTurbulenceSolve ? 2 : 1),"Lost grouped mathematical provenance.");
                operations.bind(node.operation,node.provider,[&,id=node.operation]{++counts[id];order.push_back(id);});
            });
            SF::Run::PlanExecutor::execute(resolved.solvePlan,operations);
            require(counts[OpIds::EeTurbulencePrepare]==2 && counts[OpIds::EeTurbulenceSolve]==(ras ? 2 : 0)
                && counts[OpIds::EeInterphaseCompute]==2 && counts[OpIds::EeSourcesAssemble]==2
                && counts[OpIds::EeEnergySolve]==2 && counts[OpIds::EePressureSolve]==8
                && counts[OpIds::EeStepCommit]==1,"A kernel repeated per phase or escaped its physical-step/outer contract.");
            for (std::size_t i=0;i<order.size();++i) {
                if (order[i]==OpIds::EeTurbulencePrepare) require(i>0 && order[i-1]==OpIds::EeSourcesAssemble && order[i+1]==OpIds::EeSourcesValidate,"Prepare moved relative to sources/predictor.");
                if (order[i]==OpIds::EeEnergySolve) require(order[i+1]==(ras ? OpIds::EeTurbulenceSolve : OpIds::EeBoundaryFinal),"Transport moved relative to energy/final boundary.");
            }
            auto parallel=request;parallel.parallel=true;const auto unsupported=build(config,parallel);
            require(unsupported.runtime.report.status==RuntimeStatus::Unsupported,"Serial phase transport advertised MPI capability.");
            const auto reject=[&](const ExecutionProgram& program,const ExecutableEquationSystem& equations,const std::vector<NumericalBinding>& numerics) {
                bool failed=false;try {(void)compileExecutionProgram(equations,program,numerics,builtinProviders());} catch (const std::exception&) {failed=true;}
                require(failed,"Invalid native Eulerian turbulence inherited a runnable kernel.");
            };
            auto bad=resolved.solvePlan.sourceProgram;bad.root.children[0].children.erase(bad.root.children[0].children.begin());
            reject(bad,resolved.executableSystem,resolved.numericalSelection.bindings);
            bad=resolved.solvePlan.sourceProgram;std::swap(bad.root.children[0].children[0],bad.root.children[0].children[1]);
            reject(bad,resolved.executableSystem,resolved.numericalSelection.bindings);
            if (ras) {
                bad=resolved.solvePlan.sourceProgram;bad.root.children[0].children.back().children.pop_back();
                reject(bad,resolved.executableSystem,resolved.numericalSelection.bindings);
                bad=resolved.solvePlan.sourceProgram;auto& group=bad.root.children[0].children.back();group.children.push_back(group.children.front());
                reject(bad,resolved.executableSystem,resolved.numericalSelection.bindings);
                bad=resolved.solvePlan.sourceProgram;std::swap(bad.root.children[0].children[4],bad.root.children[0].children[5]);
                reject(bad,resolved.executableSystem,resolved.numericalSelection.bindings);
                auto math=resolved.executableSystem;auto equation=math.registry.at("k."+selected[0]);
                equation.lhs=FormulaExpr::op("ddt",{FormulaExpr::symbol("k."+selected[0])});math.registry.replace(equation);
                reject(resolved.solvePlan.sourceProgram,math,resolved.numericalSelection.bindings);
                auto bindings=resolved.numericalSelection.bindings;
                for (auto& binding:bindings) if (binding.method=="EulerianTurbulenceTransport") binding.inputs.pop_back();
                reject(resolved.solvePlan.sourceProgram,resolved.executableSystem,bindings);
                bindings=resolved.numericalSelection.bindings;
                for (auto& binding:bindings) if (binding.method=="EulerianTurbulenceTransport") binding.method="TurbulenceTransport";
                reject(resolved.solvePlan.sourceProgram,resolved.executableSystem,bindings);
            }
            auto badMath=resolved.executableSystem;auto closure=badMath.registry.at("mu_t."+selected[0]);
            closure.rhs=FormulaExpr::symbol("UnsupportedModel.eddyDynamicViscosity."+selected[0]);badMath.registry.replace(closure);
            reject(resolved.solvePlan.sourceProgram,badMath,resolved.numericalSelection.bindings);
            auto wrongBacking=resolved.executableSystem;wrongBacking.state=StateRegistry{};
            for (auto symbol:resolved.executableSystem.state.symbols()) {
                if (symbol.id=="mu_t."+selected[0]) symbol.storageKey="phase99.mu_t";
                wrongBacking.state.add(symbol);
            }
            reject(resolved.solvePlan.sourceProgram,wrongBacking,resolved.numericalSelection.bindings);
            auto poisoned=resolved.solvePlan.compiledProgram;
            for (auto& call:poisoned.steps) if (eulerianTurbulenceContract(call)) {call.backendProvider="flow.eulerian-pressure";break;}
            bool integrityFailed=false;try {validateEulerianTurbulenceBindings(poisoned,request.phaseNames,tc);} catch (const std::runtime_error&) {integrityFailed=true;}
            require(integrityFailed,"Runtime accepted a changed frozen provider owner.");
            auto wrong=tc;wrong.phaseNames={"water"};if (selected==wrong.phaseNames) wrong.phaseNames={"air"};
            bool failed=false;try {validateEulerianTurbulenceBindings(resolved.solvePlan.compiledProgram,request.phaseNames,wrong);} catch (const std::runtime_error&) {failed=true;}
            require(failed,"Runtime changed phase grouping after compilation.");
            // Alias metadata is resolved by physical phase slot, not position in selected phases.
            SF::Field field;field.setup(2,2,1,1,5);SF::State::StateBundle bundle;bundle.patches={&field};
            std::vector<SF::ScalarField> backing;backing.reserve(turbulence.states.size());
            ExecutableEquationSystem stateOnly;
            for (const auto& symbol:turbulence.states) {
                stateOnly.state.add(symbol);backing.emplace_back();backing.back().setupLike(field,symbol.id,.01);
                bundle.distributed.add(SF::State::scalarView(symbol.storageKey,0,field,backing.back()));
            }
            auto realized=realizeState(stateOnly.state,{},bundle);
            for (std::size_t i=0;i<turbulence.states.size();++i) {
                const auto& symbol=turbulence.states[i];const auto slot=symbol.nameSpace=="air" ? 1 : 0;
                require(symbol.storageKey=="phase"+std::to_string(slot)+"."+symbol.id.substr(0,symbol.id.find('.')),"Selected phase index used as physical slot.");
                realized.at(symbol.id).fields.front()->write(0,0,.3+i);
                require(backing[i].values()[0]==.3+i,"STATE copied its turbulence backing.");
            }
            // Equivalent user-authored native math follows exactly the same compilation path.
            auto user=request;user.modelContributions.clear();auto authored=turbulence;authored.records.clear();
            for (auto& e:authored.registeredEquations) e.origin={OriginKind::User,"equivalent user equations"};
            for (auto& state:authored.states) state.origin={OriginKind::User,"equivalent user storage"};
            user.userContributions={authored};
            const auto equivalent=build(config,user);std::vector<std::string> a,b;
            leaves(resolved.solvePlan.root,[&](const auto& n){a.push_back(n.operation+":"+n.provider);});
            leaves(equivalent.solvePlan.root,[&](const auto& n){b.push_back(n.operation+":"+n.provider);});require(a==b,"Model/user provenance changed execution.");
        }
}
