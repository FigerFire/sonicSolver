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
    if(spec.stageConservativeTransport && (spec.eulerian || spec.model!="kOmegaSST"))throw std::runtime_error("Stage conservative turbulence currently supports single-fluid SST only.");
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
        system.placement.push_back(spec.stageConservativeTransport?System::PlacementRequirement{{},"equation:k",{"equation:energy"},{"kind:Commit"}}:System::PlacementRequirement{{},"equation:k",{},{"first-of:equation:continuity|equation:momentum"}});
        system.placement.push_back({{},"equation:"+std::string(kEpsilon ? "epsilon" : "omega"),{"equation:k"},{spec.stageConservativeTransport?"kind:Commit":"first-of:equation:continuity|equation:momentum"}});
        const std::string second=kEpsilon ? "epsilon" : "omega";
        for (const auto& id:std::vector<std::string>{"k",second,"mu_t"}) {
            System::StateSymbol symbol;
            symbol.id=id; symbol.name=id;
            symbol.role=id=="mu_t" ? System::StateRole::Derived : System::StateRole::Transported;
            symbol.storageBinding=System::StorageBinding::ProviderDistributed;
            symbol.storageKey=id;
            symbol.nameSpace="turbulence";
            if (id=="mu_t") {
                symbol.evaluation=System::StateEvaluation::Materialized;
                symbol.dependencies={"k",second,"rho","U"};
            }
            symbol.origin={System::OriginKind::Model,spec.model};
            symbol.restartEligible=false; // Manager has no restart reader for these arrays.
            symbol.initializationRequired=symbol.boundaryRequired=id!="mu_t";
            symbol.signature=mathematicalSignature(symbol.id);system.addState(std::move(symbol));
        }
        if(spec.stageConservativeTransport)for(const auto& id:{"SST.nu.k","SST.nu.omega","SST.production.k","SST.production.omega","SST.sink.k","SST.sink.omega","SST.cross","SST.F1"}) {
            System::StateSymbol coefficient;coefficient.id=coefficient.name=coefficient.storageKey=id;
            coefficient.role=System::StateRole::Derived;coefficient.outputEligible=false;coefficient.storageBinding=System::StorageBinding::TransientWorkspace;coefficient.runtimeStorageRequired=false;
            coefficient.initializationRequired=coefficient.boundaryRequired=coefficient.restartEligible=false;
            coefficient.evaluation=System::StateEvaluation::Materialized;
            if(coefficient.id=="SST.sink.k")coefficient.dependencies={"rho","omega"};
            else if(coefficient.id=="SST.sink.omega")coefficient.dependencies={"SST.F1","rho","omega"};
            else if(coefficient.id=="SST.nu.k" || coefficient.id=="SST.nu.omega") {
                coefficient.dependencies={"SST.F1","mu_t","rho"};
                coefficient.workspaceBoundary=System::WorkspaceBoundary::NearestInterior;
            }
            coefficient.signature=mathematicalSignature(coefficient.id);system.addState(std::move(coefficient));
        }
        if(spec.stageConservativeTransport) {
            if(!spec.laminarDynamicViscosity || !std::isfinite(*spec.laminarDynamicViscosity) || *spec.laminarDynamicViscosity<0)
                throw std::runtime_error("Stage SST contribution requires the configured finite nonnegative laminar dynamic viscosity.");
            const std::vector<std::pair<std::string,double>> constants={
                {"betaStar",stageSstBetaStar},{"sigmaK1",stageSstSigmaK1},{"sigmaK2",stageSstSigmaK2},
                {"sigmaW1",stageSstSigmaW1},{"sigmaW2",stageSstSigmaW2},
                {"beta1",stageSstBeta1},{"beta2",stageSstBeta2},{"laminarMu",*spec.laminarDynamicViscosity}};
            for(const auto& [name,value]:constants) {
                System::StateSymbol constant;constant.id=constant.name="kOmegaSST."+name;
                constant.signature=mathematicalSignature(constant.id);constant.constantValue=value;
                constant.runtimeStorageRequired=constant.initializationRequired=constant.boundaryRequired=constant.restartEligible=constant.outputEligible=false;
                system.addState(std::move(constant));
            }
            const std::vector<std::string> outputs={"SST.F1","mu_t","SST.production.k","SST.production.omega","SST.cross"};
            for(auto equation:stageClosureMathematics())system.addEquation(std::move(equation));
            System::ExecutionScope closure;closure.kind=System::ExecutionKind::EquationCall;
            closure.frequency=System::ExecutionFrequency::Stage;closure.step={"SST.F1",{"SST.F1"},"sst.closure"};
            system.addExecution(std::move(closure));system.bindNumerics({"SST.F1","SSTNonlocalClosure",outputs});
            system.bindInputs({"sst.closure",{{"rho",System::StateVersion::Stage},{"rhoU",System::StateVersion::Stage},{"k",System::StateVersion::Stage},{"omega",System::StateVersion::Stage}},{System::OriginKind::Model,"kOmegaSST.nonlocal"},System::ReferenceKind::Step});
            system.placement.push_back({{},"equation:SST.F1",{},{"equation:SST.sink.k","equation:SST.nu.k","equation:SST.nu.omega","equation:SST.sink.omega"}});
            std::string previous="SST.F1";
            for(auto equation:stageLocalMathematics()) {
                const auto id=equation.id,step="sst."+id.substr(4);
                system.addEquation(std::move(equation));
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
        for (const auto& id:std::vector<std::string>{"k",second}) {
            system.addEquation(spec.stageConservativeTransport?stageTransportMathematics(id):transportMathematics(second,id));
            System::ExecutionScope call;
            call.kind=System::ExecutionKind::EquationCall;
            call.step={id,{id}};
            system.addExecution(std::move(call));
            system.bindNumerics({id,spec.stageConservativeTransport?"CapacityTransportSharedFluxCentral2":"TurbulenceTransport",{"k",second}});
        }
        system.addClosure("mu_t from "+spec.model,{"k",second},{"mu_t"},{"mu_t"});
        return;
    }
    addEulerian(system,spec,true);
}

} // namespace SF::Turbulence
