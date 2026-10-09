#include "models/physics/interfaceModel/levelSet/SF_levelSetSystemContribution.h"
#include "models/physics/interfaceModel/levelSet/SF_reinit.h"
#include "solver/system/SF_systemBuilder.h"
#include "solver/run/SF_planExecutor.h"
#include "core/system/SF_operationIds.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>
using namespace SF::System;
namespace {
void require(bool condition,const char* message) {if (!condition) throw std::runtime_error(message);}
}
int main() {
    SF::FDM::SolverConfig config;
    config.numerics.timeRecipe=SF::FDM::builtInTimeRecipe(SF::FDM::TimeRecipeId::SSPRK3);
    config.numerics.recipes.time=config.numerics.timeRecipe;
    SystemContribution contribution;
    SF::Physics::InterfaceModels::LevelSetContribution::Spec spec;spec.reinitializationSteps=4;spec.pseudoTimeStep=.01;
    SF::Physics::InterfaceModels::LevelSetContribution::contribute(contribution,spec);
    require(contribution.legacyEquations.empty() && contribution.legacyExecution.empty(),"Level-set retained flat DSL/HOW authority.");
    BuildRequest request;request.composition.stateDeclared=true;request.composition.solutionVariables={"rho","rhoU","rhoE"};
    request.singleFluidPreset=SingleFluidPresetSpec{false};request.modelContributions={contribution};
    auto system=build(config,request);
    require(system.runtime.report.status==RuntimeStatus::Runnable,"Native level-set falsely unsupported.");
    SF::Run::OpRegistry operations;std::vector<std::string> trace;
    const auto bind=[&](const auto& self,const SolvePlanNode& node)->void {
        if (!node.operation.empty()) operations.bind(node.operation,node.provider,[&,op=node.operation] {trace.push_back(op);});
        for (const auto& child:node.children) self(self,child);
    };
    bind(bind,system.solvePlan.root);SF::Run::PlanExecutor::execute(system.solvePlan,operations);
    const std::vector<std::string> expected{
        OpIds::FlowStepPrepare,OpIds::FlowDtCompute,OpIds::FlowStepBegin,
        OpIds::ExplicitStageExecute,OpIds::ExplicitStageExecute,OpIds::ExplicitStageExecute,
        OpIds::FlowFinalizeBegin,OpIds::LevelSetReference,
        OpIds::LevelSetPseudoStage,OpIds::LevelSetPseudoStage,OpIds::LevelSetPseudoStage,OpIds::LevelSetPseudoStage,
        OpIds::LevelSetGeometry,OpIds::FlowStepCommit,OpIds::TimeCommit};
    require(trace==expected,"Physical RK, sign snapshot, pseudo-time, geometry or commit order changed.");
    require(std::count(trace.begin(),trace.end(),OpIds::TimeCommit)==1,"Pseudo-time advanced physical clock.");
    auto changed=contribution;changed.registeredEquations.front().rhs=FormulaExpr::constantValue(1.0);
    request.modelContributions={changed};bool rejected=false;
    try {(void)build(config,request);} catch(const std::runtime_error&) {rejected=true;}
    require(rejected,"Level-set provider silently accepted different mathematics.");

    request.parallel=true;request.modelContributions={contribution};request.capabilities.mpi=true;
    require(build(config,request).runtime.report.status!=RuntimeStatus::Runnable,"Coupled pseudo-time falsely supported.");
    request.parallel=false;spec.continuousSurfaceForce=true;spec.surfaceTension=.072;spec.interfaceThickness=.1;
    contribution={};SF::Physics::InterfaceModels::LevelSetContribution::contribute(contribution,spec);
    request.modelContributions={contribution};auto csf=build(config,request);
    require(csf.solvePlan.compiledProgram.steps[3].calls.size()==3,"CSF momentum/energy missing from fused native stage.");
    changed=contribution;
    for (auto& equation:changed.registeredEquations) if (equation.id=="interfaceEnergySource") equation.rhs=FormulaExpr::constantValue(0);
    request.modelContributions={changed};rejected=false;
    try {(void)build(config,request);} catch(const std::runtime_error&) {rejected=true;}
    require(rejected,"CSF energy exchange could drift independently of native mathematics.");

    SF::Field field;field.setup(32,1,1,4,5);
    SF::Physics::Multiphase::LevelSetField scalar;scalar.setupLike(field);
    for (int cell=0;cell<field.TotalSize();++cell) {
        int i,j,k;field.getIJK(cell,i,j,k);
        field.X(i,j,k)=.1*(i-4);field.Y(i,j,k)=j;field.Z(i,j,k)=k;
        field.XiX(i,j,k)=10;field.EtY(i,j,k)=1;field.ZeZ(i,j,k)=1;field.Jac(i,j,k)=1;
        field.CellFlag(i,j,k)=SF::FLUID_CELL;
        scalar.values()[cell]=1.4*(field.X(i,j,k)-1.55);
    }
    const auto initial=scalar.values();auto batched=scalar;
    SF::Physics::Multiphase::ReinitOptions options{12,.02,5,1e-6,2,1,{}};
    int callbacks=0;options.prepareStage=[&] {++callbacks;};
    SF::Physics::Multiphase::Reinit::advance(field,batched,options);
    require(callbacks==12,"Pseudo-time boundary/halo callback count changed.");
    SF::Physics::Multiphase::Reinit::Workspace workspace;
    SF::Physics::Multiphase::Reinit::begin(field,scalar,options,workspace);
    for (int stage=0;stage<12;++stage) SF::Physics::Multiphase::Reinit::stage(field,scalar,options,workspace);
    require(workspace.reference==initial && scalar.values()==batched.values(),"Reinitialization refroze sign or changed stage arithmetic.");
    const int ng=field.NG(),j=ng,k=ng;
    double before=0,after=0;
    for (int i=ng+8;i<ng+24;++i) {
        const auto left=field.getIdx(i-1,j,k),right=field.getIdx(i+1,j,k);
        before+=std::abs((initial[right]-initial[left])/.2-1);
        after+=std::abs((scalar.values()[right]-scalar.values()[left])/.2-1);
    }
    require(after<before,"Pseudo-time did not restore signed-distance slope.");
}
