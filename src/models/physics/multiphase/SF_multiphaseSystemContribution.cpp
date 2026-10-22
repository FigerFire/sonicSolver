#include "SF_multiphaseSystemContribution.h"
#include <stdexcept>
#include <algorithm>
#include <map>
namespace SF::Physics::Multiphase {
namespace {
using E=System::FormulaExpr;
Math::Signature signature(const std::string& name) {
    if(name=="U")return Math::signature(1,{0,1,-1,0,0,0,0});
    if(name=="thermal.T")return Math::signature(0,{0,0,0,1,0,0,0});
    if(name=="rhoU")return Math::signature(1,{1,-2,-1,0,0,0,0});
    if(name=="rhoE")return Math::signature(0,{1,-1,-2,0,0,0,0});
    if(name=="momentumFlux")return Math::signature(2,{1,-1,-2,0,0,0,0});
    if(name=="energyFlux")return Math::signature(1,{1,0,-3,0,0,0,0});
    if(name.rfind("componentFlux.",0)==0 || name.rfind("phaseMassFlux.",0)==0)return Math::signature(1,{1,-2,-1,0,0,0,0});
    if(name.rfind("partialDensity.",0)==0 || name.rfind("phaseMass.",0)==0)return Math::signature(0,{1,-3,0,0,0,0,0});
    if(name=="phaseChangeMdot" || name=="phaseMassTransfer")return Math::signature(0,{1,-3,-1,0,0,0,0});
    if(name=="phaseMassDiffusivity")return Math::signature(0,{0,2,-1,0,0,0,0});
    throw std::runtime_error("Unknown multiphase mathematical symbol: "+name);
}
E symbol(const std::string& name) {return E::symbol(name,signature(name));}

void addState(System::SystemContribution& system,const std::string& name,
        const std::string& storage,int offset,int components,bool primary) {
    System::StateSymbol state;state.id=state.name=name;state.storageKey=storage;state.componentOffset=offset;
    state.components=components;state.shape=components==1 ? System::ValueShape::Scalar : System::ValueShape::Vector;
    state.role=primary ? System::StateRole::Primary : System::StateRole::Transported;
    state.storageBinding=storage=="conservative" ? System::StorageBinding::PackedDistributed : System::StorageBinding::NamedDistributed;
    state.signature=signature(name);state.restartEligible=false;system.addState(std::move(state));
}
void addBalance(System::SystemContribution& system,const std::string& id,const std::string& unknown,
        const std::string& flux,const std::string& method,const std::vector<std::string>& members,int order,E source=E::constantValue(0.0)) {
    system.addEquation(System::Equation{id,E::add(E::op(Math::Op::Ddt,{symbol(unknown)}),E::op(Math::Op::Div,{symbol(flux)})),
        std::move(source),{}});
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
            E::multiply(E::constantValue(transferSigns.at(components[index])),E::op(Math::Op::Source,{symbol("phaseChangeMdot")})) : E::constantValue(0.0));
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
void contributeThermal(System::SystemContribution& system,const MultiPhaseConfig& config,bool ilw,int order) {
    const auto& t=config.temperature;
    if (!t.enabled) throw std::runtime_error("Thermal equation requires an enabled transported temperature.");
    system.recordContribution("model.thermal","independent transported temperature; does not replace EOS T");
    addState(system,"thermal.T",t.fieldName,0,1,false);
    if(t.transportEnabled)system.requireState("U");
    auto lhs=E::op(Math::Op::Ddt,{symbol("thermal.T")});
    if (t.transportEnabled) lhs=E::add(std::move(lhs),E::op(Math::Op::Div,{symbol("U"),symbol("thermal.T")}));
    auto rhs=t.diffusionEnabled ? E::op(Math::Op::Diffusion,{E::constantValue(t.diffusivity,Math::signature(0,{0,2,-1,0,0,0,0})),symbol("thermal.T")}) : E::constantValue(0);
    system.addEquation(System::Equation{"thermalTemperature",std::move(lhs),std::move(rhs),{}});
    System::ExecutionScope call;call.kind=System::ExecutionKind::EquationCall;call.order=20;
    call.step={"thermalTemperature",{"thermal.T"}};system.addExecution(std::move(call));
    const auto mode=normalizeModelType(t.boundMode);
    const double boundMode=(mode=="none" || mode=="off") ? 0. :
        (mode=="diagnose" || mode=="warn" || mode.empty()) ? 1. :
        (mode=="project" || mode=="clip" || mode=="bounded") ? 2. : -1.;
    if(boundMode<0)throw std::runtime_error("Unknown thermal boundMode: "+t.boundMode);
    System::NumericalBinding binding{"thermalTemperature","ThermalEulerTransport"};
    binding.parameters={{"convection",t.transportEnabled?1.:0.},{"diffusion",t.diffusionEnabled?1.:0.},
        {"D",t.diffusivity},{"bounded",t.bounded?1.:0.},{"lower",t.lowerBound},{"upper",t.upperBound},
        {"boundMode",boundMode},
        {"ilw",ilw?1.:0.},{"ilwOrder",(double)order}};
    system.bindNumerics(std::move(binding));
}
void contributeMixture(System::SystemContribution& system,const MultiPhaseConfig& config) {
    if (config.alpha.phaseName.empty()) throw std::runtime_error("Mixture contribution requires a tracked phase.");
    const auto mass="phaseMass."+config.alpha.phaseName;
    system.recordContribution("model.mixture","tracked conservative phase mass, derived alpha and thermal closure");
    system.requireProvider("equation.mixture","bind tracked phase-mass state and closure ports");
    addState(system,mass,mass,0,1,false);
    auto lhs=E::op(Math::Op::Ddt,{symbol(mass)});
    if (config.alpha.transportEnabled) lhs=E::add(std::move(lhs),E::op(Math::Op::Div,{symbol("phaseMassFlux."+config.alpha.phaseName)}));
    if (config.alpha.diffusionEnabled) lhs=E::add(std::move(lhs),E::op(Math::Op::Diffusion,{symbol("phaseMassDiffusivity"),symbol(mass)}));
    system.addEquation(System::Equation{mass,std::move(lhs),config.phaseChange.enabled ? E::op(Math::Op::Source,{symbol("phaseMassTransfer")}) : E::constantValue(0.0),{}});
    System::ExecutionScope call;call.kind=System::ExecutionKind::EquationCall;call.order=70;call.step={mass,{mass}};
    system.addExecution(std::move(call));system.bindNumerics({mass,"MixtureBalance"});
    system.addClosure("alpha = phaseMass/rho_tracked; alpha is a derived cache, never a second integrated quantity");
    if (config.temperature.enabled) system.addClosure("T is refreshed from conservative energy and phaseMass after stage/commit; no independent temperature advancement for mixture");
}
}
