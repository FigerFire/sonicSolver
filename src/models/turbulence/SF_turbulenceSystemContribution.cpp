/// @file SF_turbulenceSystemContribution.cpp
/// @brief Turbulence unknown, equation, closure, and schedule contributions.

#include "SF_turbulenceSystemContribution.h"
#include "SF_transportMathematics.h"
#include "SF_eulerianTransportMathematics.h"
#include <algorithm>
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
    system.addState(std::move(symbol));
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
    const auto contributeCalls=[&](const auto& members,const char* method,int order) {
        for (const auto& id:members) {
            System::ExecutionScope call;call.kind=System::ExecutionKind::EquationCall;
            call.order=order;call.step={id,{id}};system.addExecution(std::move(call));
            system.bindNumerics({id,method,members});
        }
    };
    contributeCalls(closures,"EulerianTurbulenceClosure",40);
    if (transport) contributeCalls(inputs,"EulerianTurbulenceTransport",50);
    system.addClosure("mu_t from "+spec.model);
}

} // namespace

void contribute(
        System::SystemContribution& system,
        const SystemContributionSpec& spec) {
    system.recordContribution(
        "model.turbulence","turbulence equation/closure contribution");
    system.requireProvider("closure.turbulence",
                           "bind turbulence closure and transported equations");
    const bool kEpsilon = spec.model == "kEpsilon";
    const bool kOmega = spec.model == "kOmegaSST";
    if (!kEpsilon && !kOmega) {
        system.addClosure(
            spec.model.empty() ? "turbulence closure" : "turbulence: "+spec.model);
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
            system.addState(std::move(closure));
            system.addEquation(System::Equation{"mu_t",System::FormulaExpr::symbol("mu_t"),
                System::FormulaExpr::symbol("Smagorinsky.eddyDynamicViscosity"),
                {System::OriginKind::Model,spec.model},true});
            System::ExecutionScope call;
            call.kind=System::ExecutionKind::EquationCall;call.order=1;
            call.step={"mu_t",{"mu_t"}};
            system.addExecution(std::move(call));
            system.bindNumerics({"mu_t","TurbulenceClosure"});
        }
        return;
    }
    if (!spec.eulerian) {
        system.requireState("rho");system.requireState("U");
        const std::string second=kEpsilon ? "epsilon" : "omega";
        for (const auto& id:std::vector<std::string>{"k",second,"mu_t"}) {
            System::StateSymbol symbol;
            symbol.id=id; symbol.name=id;
            symbol.role=id=="mu_t" ? System::StateRole::Derived : System::StateRole::Transported;
            symbol.storageBinding=System::StorageBinding::ProviderDistributed;
            symbol.storageKey=id;
            symbol.nameSpace="turbulence";
            symbol.origin={System::OriginKind::Model,spec.model};
            symbol.restartEligible=false; // Manager has no restart reader for these arrays.
            symbol.initializationRequired=symbol.boundaryRequired=id!="mu_t";
            system.addState(std::move(symbol));
        }
        for (const auto& id:std::vector<std::string>{"k",second}) {
            system.addEquation(transportMathematics(second,id));
            System::ExecutionScope call;
            call.kind=System::ExecutionKind::EquationCall;
            call.order=id=="k" ? 1 : 2;
            call.step={id,{id}};
            system.addExecution(std::move(call));
            system.bindNumerics({id,"TurbulenceTransport",{"k",second}});
        }
        system.addClosure("mu_t from "+spec.model);
        return;
    }
    addEulerian(system,spec,true);
}

} // namespace SF::Turbulence
