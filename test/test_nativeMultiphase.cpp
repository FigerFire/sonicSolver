#include "models/physics/multiphase/SF_multiphaseSystemContribution.h"
#include "solver/system/SF_systemBuilder.h"
#include <stdexcept>
using namespace SF::System;
namespace {
void require(bool value,const char* message) {if (!value) throw std::runtime_error(message);}
}
int main() {
    SF::Physics::Multiphase::MultiPhaseConfig config;
    SF::Physics::Multiphase::PhaseProperties liquid,gas;
    liquid.name="water";liquid.role=SF::Physics::Multiphase::PhaseRole::Liquid;
    gas.name="air";gas.role=SF::Physics::Multiphase::PhaseRole::Gas;
    config.phases={liquid,gas};config.phaseChange.enabled=true;
    SystemContribution contribution;SF::Physics::Multiphase::contributeHomogeneous(contribution,config);
    require(contribution.legacyEquations.empty() && contribution.legacyExecution.empty(),"Homogeneous retained flat authority.");
    require(contribution.registeredEquations.size()==4,"Homogeneous balances do not match packed components.");
    require(contribution.states[2].componentOffset==2 && contribution.states[3].componentOffset==5,"Homogeneous momentum/energy layout differs from storage.");
    require(dependsOn(contribution.registeredEquations[0].rhs,"phaseChangeMdot") && dependsOn(contribution.registeredEquations[1].rhs,"phaseChangeMdot"),"Phase-transfer sources missing in native WHAT.");
    BuildRequest request;request.homogeneousThermodynamics=true;request.modelContributions={contribution};
    SF::FDM::SolverConfig numerical;
    auto result=build(numerical,request);
    require(result.runtime.report.status!=RuntimeStatus::Runnable,"Unimplemented generic-EOS flux falsely runnable.");
    require(result.solvePlan.compiledProgram.steps.size()==4,"Homogeneous components missing from HOW.");
    config.phaseChange.enabled=false;config.alpha.phaseName="water";config.alpha.transportEnabled=true;
    contribution={};SF::Physics::Multiphase::contributeMixture(contribution,config);
    require(contribution.legacyEquations.empty() && contribution.legacyExecution.empty(),"Mixture retained flat authority.");
    require(contribution.states.size()==1 && contribution.states[0].id=="phaseMass.water","Mixture integrates alpha cache instead of phase mass.");
    request={};request.composition.stateDeclared=true;request.composition.solutionVariables={"rho","rhoU","rhoE"};
    request.singleFluidPreset=SingleFluidPresetSpec{false};request.modelContributions={contribution};
    result=build(numerical,request);
    require(result.runtime.report.status!=RuntimeStatus::Runnable,"Missing mixture thermodynamic/closure binding falsely runnable.");
    require(result.solvePlan.compiledProgram.steps.back().equationMethod=="MixtureBalance","Mixture missing native RK occurrence.");
}
