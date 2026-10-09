#include "SF_multiphaseSystemContribution.h"
#include <stdexcept>
#include <algorithm>
#include <map>
namespace SF::Physics::Multiphase {
namespace {
using E=System::FormulaExpr;
void addState(System::SystemContribution& system,const std::string& name,
        const std::string& storage,int offset,int components,bool primary) {
    System::StateSymbol state;state.id=state.name=name;state.storageKey=storage;state.componentOffset=offset;
    state.components=components;state.shape=components==1 ? System::ValueShape::Scalar : System::ValueShape::Vector;
    state.role=primary ? System::StateRole::Primary : System::StateRole::Transported;
    state.storageBinding=storage=="conservative" ? System::StorageBinding::PackedDistributed : System::StorageBinding::NamedDistributed;
    state.restartEligible=false;system.addState(std::move(state));
}
void addBalance(System::SystemContribution& system,const std::string& id,const std::string& unknown,
        const std::string& flux,const std::string& method,const std::vector<std::string>& members,int order,E source=E::constantValue(0.0)) {
    system.addEquation(System::Equation{id,E::add(E::op("ddt",{E::symbol(unknown)}),E::op("div",{E::symbol(flux)})),
        std::move(source),{},true});
    System::ExecutionScope call;call.kind=System::ExecutionKind::EquationCall;call.order=order;
    call.step={id,{unknown}};system.addExecution(std::move(call));system.bindNumerics({id,method,members});
}
}
void contributeHomogeneous(System::SystemContribution& system,const MultiPhaseConfig& config) {
    system.recordContribution("model.homogeneousMultiphase","partial-density conservative equation bundle");
    system.requireProvider("thermodynamics.homogeneous","one authoritative component EOS and packed state layout");
    std::vector<std::string> components;
    for (const auto& phase:config.phases) {
        if (phase.species.empty()) components.push_back("partialDensity."+phase.name+".bulk");
        else for (const auto& species:phase.species) components.push_back("partialDensity."+phase.name+"."+species.name);
    }
    if (components.empty()) throw std::runtime_error("Homogeneous contribution requires explicit components.");
    std::map<std::string,double> transferSigns;
    if (config.phaseChange.enabled) {
        const auto route=[&](PhaseRole role,const std::string& key,double sign) {
            for (const auto& phase:config.phases) if (phase.role==role) {
                const auto selection=config.phaseChange.selections.find(key);
                const auto species=selection==config.phaseChange.selections.end() ? std::string() : selection->second;
                if (species.empty() && phase.species.size()>1)
                    throw std::runtime_error("Homogeneous phase transfer requires an explicit species for "+phase.name);
                const auto name="partialDensity."+phase.name+"."+(species.empty() ?
                    (phase.species.empty() ? "bulk" : phase.species.front().name) : species);
                if (std::find(components.begin(),components.end(),name)==components.end())
                    throw std::runtime_error("Homogeneous transfer species is not registered: "+name);
                transferSigns[name]=sign;return;
            }
            throw std::runtime_error("Homogeneous phase transfer requires explicit liquid and gas roles.");
        };
        route(PhaseRole::Liquid,"liquidSpecies",-1.0);route(PhaseRole::Gas,"vaporSpecies",1.0);
    }
    auto members=components;members.push_back("rhoU");members.push_back("rhoE");
    for (std::size_t index=0;index<components.size();++index) {
        addState(system,components[index],"conservative",(int)index,1,true);
        addBalance(system,components[index],components[index],"componentFlux."+components[index],"HomogeneousBalance",members,10+(int)index,transferSigns.count(components[index]) ?
            E::multiply(E::constantValue(transferSigns.at(components[index])),E::op("source",{E::symbol("phaseChangeMdot")})) : E::constantValue(0.0));
    }
    addState(system,"rhoU","conservative",(int)components.size(),3,true);
    addState(system,"rhoE","conservative",(int)components.size()+3,1,true);
    addBalance(system,"momentum","rhoU","momentumFlux","HomogeneousBalance",members,1000);
    addBalance(system,"energy","rhoE","energyFlux","HomogeneousBalance",members,1001);
    system.requireState("U");system.requireState("p");system.requireState("T");
    system.addClosure("rho is the sum of partialDensity components; phase masses and alpha are derived, not independently transported");
    system.addClosure("shared p/T/U from the configured component equations of state");
    if (config.phaseChange.enabled) system.addClosure("phase-transfer ledger contributes opposite component mass sources and conservative energy mapping; native source lowering pending");
    System::ExecutionScope commit;commit.kind=System::ExecutionKind::Commit;commit.id="physicalStep.commit";commit.order=1000000;
    system.addExecution(std::move(commit));
}
void contributeMixture(System::SystemContribution& system,const MultiPhaseConfig& config) {
    if (config.alpha.phaseName.empty()) throw std::runtime_error("Mixture contribution requires a tracked phase.");
    const auto mass="phaseMass."+config.alpha.phaseName;
    system.recordContribution("model.mixture","tracked conservative phase mass, derived alpha and thermal closure");
    system.requireProvider("equation.mixture","bind tracked phase-mass state and closure ports");
    addState(system,mass,mass,0,1,false);
    auto lhs=E::op("ddt",{E::symbol(mass)});
    if (config.alpha.transportEnabled) lhs=E::add(std::move(lhs),E::op("div",{E::symbol("phaseMassFlux."+config.alpha.phaseName)}));
    if (config.alpha.diffusionEnabled) lhs=E::add(std::move(lhs),E::op("diffusion",{E::symbol("phaseMassDiffusivity"),E::symbol(mass)}));
    system.addEquation(System::Equation{mass,std::move(lhs),config.phaseChange.enabled ? E::op("source",{E::symbol("phaseMassTransfer")}) : E::constantValue(0.0),{},true});
    System::ExecutionScope call;call.kind=System::ExecutionKind::EquationCall;call.order=70;call.step={mass,{mass}};
    system.addExecution(std::move(call));system.bindNumerics({mass,"MixtureBalance"});
    system.addClosure("alpha = phaseMass/rho_tracked; alpha is a derived cache, never a second integrated quantity");
    if (config.temperature.enabled) system.addClosure("T is refreshed from conservative energy and phaseMass after stage/commit; no independent temperature advancement for mixture");
}
}
