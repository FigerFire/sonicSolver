/// @file SF_turbulenceSystemContribution.cpp
/// @brief Turbulence unknown, equation, closure, and schedule contributions.

#include "SF_turbulenceSystemContribution.h"
#include "SF_transportMathematics.h"
#include "core/system/SF_scheduleIds.h"

#include "core/system/SF_systemContribution.h"

#include <stdexcept>
#include <utility>

namespace SF::Turbulence {
namespace {

void addScalar(
        System::SystemContribution& system,
        const std::string& id,
        const std::string& name,
        const std::string& storage,
        const std::string& nameSpace) {
    System::StateSymbol unknown;
    unknown.id = id;
    unknown.name = name;
    unknown.components = 1;
    unknown.shape = System::ValueShape::Scalar;
    unknown.role = System::StateRole::Transported;
    unknown.storageBinding = System::StorageBinding::NamedDistributed;
    unknown.storageKey = storage;
    unknown.nameSpace = nameSpace;
    system.addState(std::move(unknown));
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
        if (!spec.eulerian && spec.model=="Smagorinsky") {
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
    // Legacy Eulerian compatibility: its assembly and shared-pressure policy remain unchanged.
    system.requireProvider("equation.turbulence-transport",
                           "bind Eulerian transported turbulence state and equations");

    std::vector<std::string> phases = spec.phases;
    if (!spec.eulerian) phases = {""};
    else if (phases.empty()) {
        throw std::runtime_error(
            "Eulerian transported turbulence equations require explicit phases.");
    }
    for (std::size_t phaseIndex = 0; phaseIndex < phases.size(); ++phaseIndex) {
        const std::string& phase = phases[phaseIndex];
        const std::string suffix = phase.empty() ? "" : "."+phase;
        const std::string k = "k"+suffix;
        const std::string second = (kEpsilon ? "epsilon" : "omega")+suffix;
        const std::string kEquation = "k"+suffix;
        const std::string secondEquation =
            (kEpsilon ? "epsilon" : "omega")+suffix;
        const std::string storage = phase.empty()
            ? "" : "phase"+std::to_string(phaseIndex)+".";
        addScalar(system,k,"turbulent kinetic energy"+
            (phase.empty() ? "" : " "+phase),storage+"k",
            phase.empty() ? "turbulence" : phase);
        addScalar(system,second,
            std::string(kEpsilon ? "turbulence dissipation "
                                 : "specific dissipation ")+phase,
            storage+(kEpsilon ? "epsilon" : "omega"),
            phase.empty() ? "turbulence" : phase);
        system.addEquation(
            {kEquation,"turbulent kinetic energy transport","conservation",{k}},
            SF::Equation::named(kEquation,
                SF::Equation::ddt({k}) + SF::Equation::div({"flux."+k})
                    + SF::Equation::diffusion({"diffusivity."+k},{k})
                    == SF::Equation::Symbol{"source."+k}));
        system.addEquation(
            {secondEquation,kEpsilon ? "epsilon transport" : "omega transport",
             "conservation",{second}},
            SF::Equation::named(secondEquation,
                SF::Equation::ddt({second}) + SF::Equation::div({"flux."+second})
                    + SF::Equation::diffusion(
                        {"diffusivity."+second},{second})
                    == SF::Equation::Symbol{"source."+second}));
        for (const auto& item:std::vector<std::pair<std::string,int>>{{kEquation,50},{secondEquation,51}}) {
            System::ExecutionScope call;
            call.kind=System::ExecutionKind::EquationCall;
            call.order=item.second;
            call.step={item.first,{item.first}};
            system.addLegacyExecution(std::move(call));
        }

    }
    system.addClosure("mu_t from "+spec.model);
}

} // namespace SF::Turbulence
