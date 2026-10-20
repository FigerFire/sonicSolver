#include "models/physics/multiphase/SF_multiphaseSystemContribution.h"
#include "solver/system/SF_systemBuilder.h"
#include "solver/system/SF_builtinState.h"
#include "models/physics/fluidStateModel/SF_factory.h"
#include "solver/system/SF_systemValidator.h"
#include "solver/algorithm/SF_singleFluidStepper.h"
#include "methods/numerics/scalar/SF_scalarTransport.h"
#include <iostream>
#include <stdexcept>
using namespace SF;
void require(bool c,const char* s){if(!c)throw std::runtime_error(s);}
int main(){
 Math::deactivateDirection(2);
 for(bool ilw:{false,true})for(bool convection:{false,true})for(bool diffusion:{false,true}) {
    FDM::SolverConfig config;config.numerics.timeRecipe=FDM::builtInTimeRecipe(FDM::TimeRecipeId::ForwardEuler);
    config.numerics.maxDeltaT=.001;
    Physics::Multiphase::MultiPhaseConfig settings;settings.enabled=true;settings.type="thermal";
    auto& t=settings.temperature;t.enabled=true;t.fieldName="independentTemperature";t.bounded=false;
    t.transportEnabled=convection;t.diffusionEnabled=diffusion;t.diffusivity=.2;
    System::SystemContribution contribution;Physics::Multiphase::contributeThermal(contribution,settings,ilw,ilw?3:0);
    if(convection) {
        contribution.addState(System::builtinState(System::BuiltinState::Velocity));
        for(auto variable:{System::BuiltinState::Density,System::BuiltinState::Momentum}) {
            auto symbol=System::builtinState(variable);symbol.storageBinding=System::StorageBinding::PackedDistributed;
            symbol.runtimeStorageRequired=true;contribution.addState(symbol);
        }
    }
    System::ExecutionScope commit;commit.kind=System::ExecutionKind::Commit;commit.id="commit";commit.order=1000000;contribution.addExecution(commit);
    System::BuildRequest request;request.includeDefaultFluidPreset=false;request.composition.stateDeclared=true;
    request.composition.solutionVariables={"thermal.T"};request.userContributions={contribution};
    auto system=System::build(config,request);System::validate(system);
    require(system.runtime.report.status==System::RuntimeStatus::Runnable,"Thermal plan not runnable");
    Field f;f.setup(9,9,1,3,5);ScalarField actual,reference;actual.setupLike(f,t.fieldName);
    for(int k=0;k<f.MZ();++k)for(int j=0;j<f.MY();++j)for(int i=0;i<f.MX();++i){
      f.X(i,j,k)=.125*(i-3);f.Y(i,j,k)=j-3;f.Z(i,j,k)=k-3;
      f.XiX(i,j,k)=8;f.EtY(i,j,k)=f.ZeZ(i,j,k)=1;f.Jac(i,j,k)=1;
      f.CellFlag(i,j,k)=FLUID_CELL;f(i,j,k,0)=1;f(i,j,k,1)=.25;f(i,j,k,4)=2.53125;
      actual(i,j,k)=300+3*std::sin(f.X(i,j,k));
    }
    std::vector<int> left,right;
    for(int j=3;j<12;++j){left.push_back(f.getIdx(3,j,3));right.push_back(f.getIdx(11,j,3));}
    f.setBoundarySets({{"Left",left},{"Right",right}});
    std::vector<BCSetting<double>> boundaries={{"Left",FIXED_VALUE,301},{"Right",ZERO_GRADIENT,0}};
    f.setStateModel(Physics::FluidStateModel::makeSingleFluidPerfectGas(1.4,1,0,.72));
    reference=actual;
    FDM::Scalar::TransportConfig numerical;numerical.name=t.fieldName;numerical.convectionEnabled=convection;
    numerical.diffusionEnabled=diffusion;numerical.diffusivity=.2;
    numerical.ilwBoundaryEnabled=ilw;numerical.ilwAccuracyOrder=ilw?3:0;
    State::StateBundle bundle;bundle.patches={&f};bundle.stateModel=f.stateModel();bundle.time=.3;bundle.step=17;
    State::VariableDescriptor variable;variable.name=t.fieldName;variable.updatePolicy=State::UpdatePolicy::ExternalExplicit;
    bundle.transported.add(variable,actual);bundle.distributed.add(State::scalarView(t.fieldName,0,f,actual));
    if(convection)bundle.distributed.add(State::conservativeView("conservative",0,f));
    FDM::ThermalEquationData data{&actual,boundaries};FDM::SolverServices services;services.thermalEquation=&data;
    SolverAlgorithm::SingleFluidStepper solver(config,system.executableSystem,system.numericalSystem,system.solvePlan,system.runtime);
    solver.bindServices(services);solver.bindSolvePlan(system.solvePlan);
    FDM::SolverState state;state.bundle=&bundle;state.maximumTimeStep=.001;solver.prepare(state);
    for(int step=0;step<5;++step){
      FDM::Scalar::advanceForwardEuler(f,reference,boundaries,numerical,.001);
      auto result=solver.advance(state);require(result.accepted && result.dt==.001,"Thermal step dt/acceptance");
      require(actual.values()==reference.values(),"First thermal difference from original Forward Euler kernel");
      require(bundle.step==18+step && std::abs(bundle.time-(.301+.001*step))<1e-14,"Thermal double update/clock/restart");
    }
    auto invalid=contribution;invalid.registeredEquations[0].rhs=System::FormulaExpr::constantValue(123);
    request.userContributions={invalid};bool rejected=false;try{(void)System::build(config,request);}catch(const std::exception&){rejected=true;}
    require(rejected,"Thermal recipe silently accepts modified WHAT");
    if(!ilw && !convection && diffusion) {
        auto bad=contribution;bad.states.clear();request.userContributions={bad};
        rejected=false;try{(void)System::build(config,request);}catch(const std::exception&){rejected=true;}
        require(rejected,"Thermal missing STATE accepted");
        bad=contribution;bad.execution.insert(bad.execution.begin(),bad.execution.front());request.userContributions={bad};
        rejected=false;try{(void)System::build(config,request);}catch(const std::exception&){rejected=true;}
        require(rejected,"Thermal duplicate HOW target accepted");
        bad=contribution;bad.numerics.front().method="UnknownThermalMethod";request.userContributions={bad};
        rejected=false;try{(void)System::build(config,request);}catch(const std::exception&){rejected=true;}
        require(rejected,"Thermal unknown WHICH silently substituted");
        request.userContributions={contribution};request.parallel=true;request.capabilities.mpi=true;
        rejected=false;try{(void)System::build(config,request);}catch(const std::exception&){rejected=true;}
        require(rejected,"Thermal MPI falsely Runnable");request.parallel=false;
    }
    request.userContributions={contribution};config.numerics.timeRecipe=FDM::builtInTimeRecipe(FDM::TimeRecipeId::ClassicalRK4);
    rejected=false;try{(void)System::build(config,request);}catch(const std::exception&){rejected=true;}
    require(rejected,"Thermal falsely supports RK4");
 }
 std::cout<<"8 thermal operators with fixedValue/zeroGradient algebraic or ILW3 BC x 5 physical steps: exact old-kernel equality; nonzero clock/step and negative WHAT/RK4 passed\n";
}
