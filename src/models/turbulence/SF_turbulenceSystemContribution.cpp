/// @file SF_turbulenceSystemContribution.cpp
/// @brief Turbulence unknown, equation, closure, and schedule contributions.

#include "SF_turbulenceSystemContribution.h"
#include "SF_transportMathematics.h"
#include "SF_eulerianTransportMathematics.h"
#include <algorithm>
#include <cmath>
#include "core/system/SF_scheduleIds.h"

#include "core/system/SF_systemContribution.h"

#include <stdexcept>
#include <utility>

namespace SF::Turbulence {
namespace {

void addEulerianState(System::SystemContribution& system,const std::string& id,
        const std::string& storage,const std::string& phase,bool closure) {
    System::StateSymbol symbol;
    symbol.id=symbol.name=id;symbol.storageKey=storage;symbol.nameSpace=phase;
    symbol.role=closure ? System::StateRole::Derived : System::StateRole::Transported;
    symbol.storageBinding=System::StorageBinding::ProviderDistributed;
    symbol.initializationRequired=symbol.boundaryRequired=!closure;
    symbol.restartEligible=false;
    symbol.signature=mathematicalSignature(symbol.id);system.addState(std::move(symbol));
}
void addEulerian(System::SystemContribution& system,const SystemContributionSpec& spec,bool transport) {
    if (spec.phases.empty() || spec.phaseOrder.empty())
        throw std::runtime_error("Eulerian turbulence requires selected phases and complete phase storage order.");
    std::vector<std::string> inputs,closures;
    const auto second=spec.model=="kEpsilon" ? "epsilon" : "omega";
    for (std::size_t i=0;i<spec.phases.size();++i) {
        const auto& phase=spec.phases[i];
        const auto at=std::find(spec.phaseOrder.begin(),spec.phaseOrder.end(),phase);
        if (phase.empty() || at==spec.phaseOrder.end()
            || std::find(spec.phases.begin(),spec.phases.begin()+i,phase)!=spec.phases.begin()+i)
            throw std::runtime_error("Invalid/duplicate selected Eulerian turbulence phase: "+phase);
        const auto prefix="phase"+std::to_string(at-spec.phaseOrder.begin())+".";
        const auto closure="mu_t."+phase;closures.push_back(closure);
        addEulerianState(system,closure,prefix+"mu_t",phase,true);
        system.addEquation(eulerianClosureMathematics(spec.model,phase));
        if (transport) for (const auto& variable:std::vector<std::string>{"k",second}) {
            const auto id=variable+"."+phase;inputs.push_back(id);
            addEulerianState(system,id,prefix+variable,phase,false);
            system.addEquation(eulerianTransportMathematics(second,phase,variable));
        }
    }
    const auto contributeCalls=[&](const auto& members,const char* method) {
        for (const auto& id:members) {
            System::ExecutionScope call;call.kind=System::ExecutionKind::EquationCall;
            call.step={id,{id}};system.addExecution(std::move(call));
            system.bindNumerics({id,method,members});
        }
    };
    contributeCalls(closures,"EulerianTurbulenceClosure");
    if (transport) contributeCalls(inputs,"EulerianTurbulenceTransport");
    system.addClosure("mu_t from "+spec.model,closures,{},closures);
}

} // namespace

void contribute(
        System::SystemContribution& system,
        const SystemContributionSpec& spec) {
    system.recordContribution(
        "model.turbulence","turbulence equation/closure contribution");
    system.requireProvider("closure.turbulence",
                           "bind turbulence closure and transported equations");
    const bool ordinary=spec.stageConservativeTransport || spec.fixedVolumeTransport;
    if(ordinary && (spec.eulerian || spec.model!="kOmegaSST"))throw std::runtime_error("Ordinary SST declarations currently support single-fluid SST only.");
    if(spec.closureDensity.empty()!=spec.closureVelocity.empty() || (!ordinary && !spec.closureVelocity.empty()))
        throw std::runtime_error("Explicit SST density/velocity inputs require paired names and the ordinary SST declaration");
    if(!ordinary && (!spec.viscosityOutput.empty() || !spec.viscosityEquation.empty()))
        throw std::runtime_error("viscosityOutput and viscosityEquation require ordinary SST declarations");
    if(!spec.viscosityEquation.empty() && !spec.fixedVolumeTransport)
        throw std::runtime_error("A replacement viscosityEquation requires fixedVolume Step: the conservative Stage residual reads the model array before final publication");
    if(spec.stressTarget && !spec.fixedVolumeTransport)
        throw std::runtime_error("A momentum stress target requires fixedVolume transport");
    const bool kEpsilon = spec.model == "kEpsilon";
    const bool kOmega = spec.model == "kOmegaSST";
    if (!kEpsilon && !kOmega) {
        if (!spec.eulerian) system.addClosure(
            spec.model.empty() ? "turbulence closure" : "turbulence: "+spec.model,{},spec.eulerian ? std::vector<std::string>{} : std::vector<std::string>{"mu_t"},spec.model=="Smagorinsky" ? std::vector<std::string>{"mu_t"} : std::vector<std::string>{});
        if (spec.eulerian && spec.model=="Smagorinsky") {
            addEulerian(system,spec,false);
        } else if (!spec.eulerian && spec.model=="Smagorinsky") {
            system.requireState("rho");system.requireState("U");
            System::StateSymbol closure;
            closure.id=closure.name=closure.storageKey="mu_t";
            closure.role=System::StateRole::Derived;
            closure.storageBinding=System::StorageBinding::ProviderDistributed;
            closure.initializationRequired=closure.boundaryRequired=closure.restartEligible=false;
            closure.origin={System::OriginKind::Model,spec.model};
            closure.signature=mathematicalSignature(closure.id);system.addState(std::move(closure));
            system.addEquation(System::Equation{"mu_t",mathematicalSymbol("mu_t"),
                mathematicalSymbol("Smagorinsky.eddyDynamicViscosity"),
                {System::OriginKind::Model,spec.model}});
            System::ExecutionScope call;
            call.kind=System::ExecutionKind::EquationCall;
            call.step={"mu_t",{"mu_t"}};
            system.addExecution(std::move(call));
            system.bindNumerics({"mu_t","TurbulenceClosure"});
            system.placement.push_back({{},"equation:mu_t",{}, {"first-of:equation:continuity|equation:momentum"}});
        }
        return;
    }
    if (!spec.eulerian) {
        system.requireState("rho");system.requireState("U");
        if(!spec.fixedVolumeTransport) {
            system.placement.push_back(spec.stageConservativeTransport?System::PlacementRequirement{{},"equation:k",{"equation:energy"},{"kind:Commit"}}:System::PlacementRequirement{{},"equation:k",{},{"first-of:equation:continuity|equation:momentum"}});
            system.placement.push_back({{},"equation:"+std::string(kEpsilon ? "epsilon" : "omega"),{"equation:k"},{spec.stageConservativeTransport?"kind:Commit":"first-of:equation:continuity|equation:momentum"}});
        }
        const std::string second=kEpsilon ? "epsilon" : "omega";
        const std::string viscosityId=spec.viscosityOutput.empty()?"mu_t":spec.viscosityOutput;
        const auto viscosity=System::FormulaExpr::symbol(viscosityId,mathematicalSignature("mu_t"));
        const std::string viscosityEquation=spec.viscosityEquation.empty()?viscosityId:spec.viscosityEquation;
        const bool builtinViscosity=spec.viscosityEquation.empty();
        const auto closureDensity=spec.closureDensity.empty()?std::string("rho"):spec.closureDensity;
        const auto closureVelocity=spec.fixedVolumeTransport && spec.closureVelocity.empty()?std::string("U"):spec.closureVelocity;
        auto localEquations=ordinary?stageLocalMathematics(viscosity):std::vector<System::Equation>{};
        for (const auto& key:std::vector<std::string>{"k",second,"mu_t"}) {
            System::StateSymbol symbol;
            symbol.id=key=="mu_t"?viscosityId:key;
            symbol.name=symbol.id;
            symbol.role=key=="mu_t" ? System::StateRole::Derived : System::StateRole::Transported;
            symbol.storageBinding=System::StorageBinding::ProviderDistributed;
            symbol.storageKey=key;
            symbol.nameSpace="turbulence";
            if (key=="mu_t") {
                symbol.evaluation=System::StateEvaluation::Materialized;
                if(builtinViscosity)
                    symbol.dependencies={"k",second,closureDensity,closureVelocity.empty()?"U":closureVelocity};
            }
            symbol.origin={System::OriginKind::Model,spec.model};
            symbol.restartEligible=false; // Manager has no restart reader for these arrays.
            symbol.initializationRequired=symbol.boundaryRequired=key!="mu_t";
            symbol.signature=mathematicalSignature(key);system.addState(std::move(symbol));
        }
        if(ordinary)for(const auto& id:{"SST.nu.k","SST.nu.omega","SST.production.k","SST.production.omega","SST.sink.k","SST.sink.omega","SST.cross","SST.F1"}) {
            System::StateSymbol coefficient;coefficient.id=coefficient.name=coefficient.storageKey=id;
            coefficient.role=System::StateRole::Derived;coefficient.outputEligible=false;coefficient.storageBinding=System::StorageBinding::TransientWorkspace;coefficient.runtimeStorageRequired=false;
            coefficient.initializationRequired=coefficient.boundaryRequired=coefficient.restartEligible=false;
            coefficient.evaluation=System::StateEvaluation::Materialized;
            if(coefficient.id=="SST.sink.k")coefficient.dependencies={"rho","omega"};
            else if(coefficient.id=="SST.sink.omega")coefficient.dependencies={"SST.F1","rho","omega"};
            else if(coefficient.id=="SST.nu.k" || coefficient.id=="SST.nu.omega") {
                coefficient.dependencies={"SST.F1",viscosityId,"rho"};
                coefficient.workspaceBoundary=System::WorkspaceBoundary::NearestInterior;
            }
            coefficient.signature=mathematicalSignature(coefficient.id);system.addState(std::move(coefficient));
        }
        if(ordinary) {
            if(!spec.laminarDynamicViscosity || !std::isfinite(*spec.laminarDynamicViscosity) || *spec.laminarDynamicViscosity<0)
                throw std::runtime_error("Stage SST contribution requires the configured finite nonnegative laminar dynamic viscosity.");
            const std::vector<std::pair<std::string,double>> constants={
                {"betaStar",stageSstBetaStar},{"sigmaK1",stageSstSigmaK1},{"sigmaK2",stageSstSigmaK2},
                {"sigmaW1",stageSstSigmaW1},{"sigmaW2",stageSstSigmaW2},
                {"omegaFloor",spec.omegaFloor},{"gamma1",5.0/9.0},{"gamma2",.44},{"a1",.31},{"beta1",stageSstBeta1},{"beta2",stageSstBeta2},{"laminarMu",*spec.laminarDynamicViscosity}};
            for(const auto& [name,value]:constants) {
                System::StateSymbol constant;constant.id=constant.name="kOmegaSST."+name;
                constant.signature=mathematicalSignature(constant.id);constant.constantValue=value;
                constant.runtimeStorageRequired=constant.initializationRequired=constant.boundaryRequired=constant.restartEligible=constant.outputEligible=false;
                system.addState(std::move(constant));
            }
            const std::vector<std::string> outputs={"SST.F1","SST.production.k","SST.production.omega","SST.cross"};
            for(auto equation:stageClosureMathematics(closureDensity,closureVelocity,viscosity,builtinViscosity))
                system.addEquation(std::move(equation));
            for(const auto& equation:localEquations)system.addEquation(equation);
            if(!spec.fixedVolumeTransport) {
                System::ExecutionScope closure;closure.kind=System::ExecutionKind::EquationCall;
                closure.frequency=System::ExecutionFrequency::Stage;closure.step={"SST.F1",{"SST.F1"},"sst.closure"};
                System::ExecutionScope viscosityCall;viscosityCall.kind=System::ExecutionKind::EquationCall;
                viscosityCall.frequency=System::ExecutionFrequency::Stage;
                viscosityCall.step={viscosityEquation,{viscosityId},"sst.viscosity"};
                system.addExecution(std::move(viscosityCall));
                system.addExecution(std::move(closure));system.bindNumerics({"SST.F1","SSTCentral2",outputs});
                std::vector<System::StateUse> reads={
                    {closureDensity,spec.closureDensity.empty()?System::StateVersion::Stage:System::StateVersion::Current},
                    {closureVelocity.empty()?"rhoU":closureVelocity,closureVelocity.empty()?System::StateVersion::Stage:System::StateVersion::Current},
                    {"k",System::StateVersion::Stage},{"omega",System::StateVersion::Stage}};
                if(builtinViscosity) {
                    system.bindNumerics({viscosityEquation,"SSTCentral2"});
                    system.bindInputs({"sst.viscosity",reads,{System::OriginKind::Model,"SST viscosity"},System::ReferenceKind::Step});
                }
                system.placement.push_back({{},"sst.viscosity",{},{"equation:SST.F1"}});
                reads.push_back({viscosityId,System::StateVersion::Stage});
                system.bindInputs({"sst.closure",std::move(reads),{System::OriginKind::Model,"kOmegaSST.nonlocal"},System::ReferenceKind::Step});
                system.placement.push_back({{},"equation:SST.F1",{},{"equation:SST.sink.k","equation:SST.nu.k","equation:SST.nu.omega","equation:SST.sink.omega"}});
                std::string previous="SST.F1";
                for(const auto& equation:localEquations) {
                    const auto id=equation.id,step="sst."+id.substr(4);
                    System::ExecutionScope coefficient;coefficient.kind=System::ExecutionKind::EquationCall;
                    coefficient.frequency=System::ExecutionFrequency::Stage;coefficient.step={id,{id},step};
                    system.addExecution(std::move(coefficient));
                    std::vector<System::StateUse> reads;
                    const auto& output=*std::find_if(system.states.begin(),system.states.end(),[&](const auto& state){return state.id==id;});
                    for(const auto& input:output.dependencies)reads.push_back({input,System::StateVersion::Stage});
                    system.bindInputs({step,std::move(reads),{System::OriginKind::Model,"kOmegaSST.stageConservative"},System::ReferenceKind::Step});
                    system.bindNumerics({id,"LocalAffine"});
                    system.placement.push_back({{},"equation:"+id,{"equation:"+previous},{"first-of:equation:continuity|equation:momentum"}});
                    previous=id;
                }
            }
        }
        for (const auto& id:std::vector<std::string>{"k",second}) {
            system.addEquation(spec.fixedVolumeTransport?fixedVolumeMathematics(id,spec.closureVelocity.empty()?"U":spec.closureVelocity):spec.stageConservativeTransport?stageTransportMathematics(id):transportMathematics(second,id));
            if(!spec.fixedVolumeTransport) {
                System::ExecutionScope call;
                call.kind=System::ExecutionKind::EquationCall;
                call.step={id,{id}};
                system.addExecution(std::move(call));
                system.bindNumerics({id,spec.stageConservativeTransport?"CapacityTransportSharedFluxCentral2":"TurbulenceTransport",{"k",second}});
            }
        }
        if(spec.fixedVolumeTransport) {
            if(spec.sharedVolumeFlux.empty())throw std::runtime_error("fixedVolume SST needs explicit sharedVolumeFlux STATE");
            const std::string velocity=spec.closureVelocity.empty()?"U":spec.closureVelocity;
            System::ExecutionScope group;group.id="sst.iteration";group.frequency=System::ExecutionFrequency::Step;
            const auto call=[&](const std::string& equation,const std::string& occurrence) {
                System::ExecutionScope node;node.kind=System::ExecutionKind::EquationCall;
                node.id=occurrence;node.step={equation,{equation},occurrence};return node;
            };
            auto viscosityCall=call(viscosityEquation,"sst.viscosity");
            viscosityCall.step.target.symbol=viscosityId;
            group.children.push_back(std::move(viscosityCall));
            if(builtinViscosity) {
                System::NumericalBinding method{viscosityEquation,"SSTCentral2"};
                method.parameters={{"constantDensity",1}};
                system.bindNumerics(std::move(method));
                system.bindInputs({"sst.viscosity",{{"rho",System::StateVersion::Current},{velocity,System::StateVersion::Current},{"k",System::StateVersion::Current},{"omega",System::StateVersion::Current}},{System::OriginKind::Model,"SST viscosity"},System::ReferenceKind::Step});
            }
            group.children.push_back(call("SST.F1","sst.closure"));
            System::NumericalBinding closureMethod{"SST.F1","SSTCentral2",{"SST.F1","SST.production.k","SST.production.omega","SST.cross"}};
            closureMethod.parameters={{"constantDensity",1}};system.bindNumerics(std::move(closureMethod));
            system.bindInputs({"sst.closure",{{"rho",System::StateVersion::Current},{velocity,System::StateVersion::Current},{"k",System::StateVersion::Current},{"omega",System::StateVersion::Current},{viscosityId,System::StateVersion::Current}},{System::OriginKind::Model,"kOmegaSST.fixedVolume"},System::ReferenceKind::Step});
            for(const auto& equation:localEquations) {
                const auto& id=equation.id;
                const auto occurrence="sst."+id.substr(4);
                group.children.push_back(call(id,occurrence));
                system.bindNumerics({id,"LocalAffine"});
                std::vector<System::StateUse> reads;
                const auto& output=*std::find_if(system.states.begin(),system.states.end(),[&](const auto& value){return value.id==id;});
                for(const auto& input:output.dependencies)reads.push_back({input,System::StateVersion::Current});
                system.bindInputs({occurrence,std::move(reads),{System::OriginKind::Model,"kOmegaSST.fixedVolume"},System::ReferenceKind::Step});
            }
            if(spec.stressTarget) {
                System::StateSymbol effective;effective.id=effective.name=effective.storageKey="SST.muEff";
                effective.signature=mathematicalSignature("mu");effective.role=System::StateRole::Derived;
                effective.storageBinding=System::StorageBinding::TransientWorkspace;effective.evaluation=System::StateEvaluation::Materialized;
                effective.dependencies={viscosityId};effective.workspaceBoundary=System::WorkspaceBoundary::NearestInterior;
                effective.runtimeStorageRequired=effective.initializationRequired=effective.boundaryRequired=effective.outputEligible=effective.restartEligible=false;
                system.addState(effective);
                using E=System::FormulaExpr;
                system.addEquation({effective.id,effective.expression(),E::add(mathematicalSymbol("kOmegaSST.laminarMu"),viscosity),{System::OriginKind::Model,"SST effective dynamic viscosity"}});
                group.children.push_back(call(effective.id,"sst.effectiveViscosity"));system.bindNumerics({effective.id,"LocalAffine"});
                system.bindInputs({"sst.effectiveViscosity",{{viscosityId,System::StateVersion::Current}},{System::OriginKind::Model,"SST stress"},System::ReferenceKind::Step});
                const auto& [equation,occurrence]=*spec.stressTarget;
                const auto stress=E::divide(E::op(Math::Op::Div,{E::multiply(effective.expression(),
                    E::op(Math::Op::DevSymmGrad,{E::symbol(velocity,mathematicalSignature("U"))}))},occurrence),mathematicalSymbol("rho"));
                system.mathematicalExtensions.push_back({{equation},stress,true,occurrence,Math::Op::Diffusion});
                System::NumericalBinding spatial;
                spatial.equation=equation;
                spatial.operators.emplace_back(equation,occurrence,"stress.central2.arithmetic-faces");
                system.bindNumerics(std::move(spatial));
            }
            System::ExecutionScope transport;transport.id="sst.transport";transport.frequency=System::ExecutionFrequency::Step;
            for(const auto& target:std::vector<std::string>{"k","omega"}) {
                const auto occurrence="sst.iterate."+target,residual="sst.defect."+target;
                System::ExecutionScope solve;solve.kind=System::ExecutionKind::While;solve.id="sst.solve."+target;
                solve.condition="res("+residual+")";solve.repetitions=100;solve.minimumIterations=1;
                solve.children.push_back(call(target,occurrence));transport.children.push_back(std::move(solve));
                System::ResidualDefinition definition;definition.id=residual;definition.occurrence=occurrence;
                definition.target=target;definition.kind=System::ResidualKind::DiscreteEquation;
                definition.origin={System::OriginKind::Model,"kOmegaSST.fixedVolume/default residual"};
                system.residuals.push_back(std::move(definition));
                System::NumericalBinding method{target,"FixedPointVolumeFluxUpwind1",{spec.sharedVolumeFlux}};
                method.operators.emplace_back(target,"diffusion."+target,"scalar.central2.arithmetic-faces");
                system.bindNumerics(std::move(method));
                std::vector<System::StateUse> reads={{target,System::StateVersion::Current},{velocity,System::StateVersion::Current},{spec.sharedVolumeFlux,System::StateVersion::Current}};
                for(const auto& input:std::vector<std::string>{"SST.nu."+target,"SST.production."+target,"SST.sink."+target})
                    reads.push_back({input,System::StateVersion::Frozen,false,"",true});
                if(target=="omega")reads.push_back({"SST.cross",System::StateVersion::Frozen,false,"",true});
                system.bindInputs({occurrence,std::move(reads),{System::OriginKind::Model,"kOmegaSST.fixedVolume"},System::ReferenceKind::Step});
            }
            if(spec.stressTarget) {
                system.addExecution(std::move(group));system.addExecution(std::move(transport));
                system.placement.push_back({"slot:coupling","sst.iteration",{},{"slot:predictor"}});
                system.placement.push_back({"slot:coupling","sst.transport",{"slot:transport"},{"first-of:slot:convergence|kind:Commit"}});
            } else {
                group.children.insert(group.children.end(),transport.children.begin(),transport.children.end());
                system.addExecution(std::move(group));
                system.placement.push_back({"slot:coupling","sst.iteration",{"slot:transport"},{"first-of:slot:convergence|kind:Commit"}});
            }
        }
        system.addClosure("mu_t from "+spec.model,{"k",second},{viscosityId},{viscosityId});
        return;
    }
    addEulerian(system,spec,true);
}

} // namespace SF::Turbulence
