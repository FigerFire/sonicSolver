#include "test/typedSymbols.h"
#include "fluid_fixture.h"
#include "solver/discretization/diffusion/SF_formulaCentral2.h"
#include "solver/algorithm/SF_flowOperations.h"
#include "solver/algorithm/scalar/SF_scalarBinding.h"
using namespace SF::SolverAlgorithm;
namespace {
System::BuildRequest sharedRequest(bool scalar=true,bool reverse=false) {
    auto r=fluidRequest(scalar,false,reverse);
    if(scalar) {
        auto& c=r.userContributions[0];c.numerics[0].method="ScalarTransportSharedFluxCentral2";
        using E=System::FormulaExpr;
        c.registeredEquations[0].rhs=E::op("diffusion",{E::multiply(SF::Test::symbol("rho"),E::constantValue(0,Math::signature(0,{0,2,-1,0,0,0,0}))),E::divide(SF::Test::symbol("rhoC"),SF::Test::symbol("rho"))});
    }
    return r;
}
void shared(FDM::TimeRecipeId id,bool constant,int nodes=41) {
    const int g=3;const double h=1./(nodes-1),dt=.0005,end=.01;
    FDM::SolverConfig config;config.numerics.maxDeltaT=dt;config.numerics.cfl=2;config.numerics.timeRecipe=FDM::builtInTimeRecipe(id);
    auto system=System::build(config,sharedRequest());System::validate(system);
    Field mesh;mesh.setup(nodes,1,1,g,5);auto model=Physics::FluidStateModel::makeSingleFluidPerfectGas(1.4,287.05,0,.72);mesh.setStateModel(model);
    ScalarField scalar;scalar.setupLike(mesh,"rhoC");
    for(int cell=0;cell<mesh.TotalSize();++cell) {
        int i,j,k;mesh.getIJK(cell,i,j,k);double x=(i-g)*h,rho=1+.1*std::sin(2*std::acos(-1.)*x);
        mesh.X(i,j,k)=x;mesh.Jac(i,j,k)=1/h;mesh.XiX(i,j,k)=1/h;mesh.EtY(i,j,k)=1;mesh.ZeZ(i,j,k)=1;
        mesh(i,j,k,0)=rho;mesh(i,j,k,1)=rho*.3;mesh(i,j,k,4)=2.5+.5*rho*.3*.3;
        scalar(i,j,k)=rho*(constant?.7:.5+.2*std::exp(-100*(x-.5)*(x-.5)));
    }
    State::StateBundle state;state.patches={&mesh};state.stateModel=model;state.registerConservativeState();state.distributed.add(State::scalarView("rhoC",0,mesh,scalar));
    auto data=std::make_shared<FDM::ScalarEquationData>();data->target="rhoC";for(auto& bc:data->boundary)bc.kind=FDM::ScalarBoundaryKind::ZeroGradient;
    Execution::LocalRuntime runtime;HomogeneousExterior exterior;FDM::SolverServices services;services.executionRuntime=&runtime;services.boundaryPipeline=&exterior;services.scalarInstances={{"advance.rhoC",data}};
    auto realized=System::realizeState(system.executableSystem.state,system.runtime,state,system.solvePlan.compiledProgram.stateViews);
    const double limit=dt;
    auto flow=std::make_shared<FlowOperations>(config,system.executableSystem,system.numericalSystem,system.solvePlan,system.runtime,state,services,realized,limit);flow->prepare();
    auto f=flow->temporalParticipant();auto scalars=scalarParticipants(system.numericalSystem,system.solvePlan,state,services,realized,limit,[flow](auto stamp){return flow->stageMassFlux(stamp);});auto& c=scalars.at(0);
    double accumulated=0,maxBudget=0,maxC=0,maxFluxDifference=0;int fluxReads=0,comparedFaces=0;
    const auto integral=[&]{double sum=0;for(int i=g;i<g+nodes;++i)sum+=scalar(i,g,g)*h;return sum;};
    double initial=integral();std::ofstream log("stage-budget-"+std::string(FDM::toString(id))+(constant?"-constant":"-pulse")+".csv");log.precision(17);log<<"step,stage,time,dt,netScalarFlux\n";
    for(int step=0;step<int(end/dt);++step) {
        state.dt=dt;f.snapshot();c.snapshot();double net=0;
        for(int stage=0;stage<config.numerics.timeRecipe.stageCount();++stage) {
            Run::ExecutionContext context;context.stageIndex=stage;context.stageCount=config.numerics.timeRecipe.stageCount();
            double t=state.time+system.numericalSystem.time.recipe.stage(stage).abscissa*dt;
            realized.beginStage(stage,context.stageCount,t);f.prepare(context);c.prepare(context);realized.stageReady();f.rhs(context);
            auto flux=flow->stageMassFlux({state.step,stage,t});double left=flux.at(0,g-1,g,g),right=flux.at(0,g+nodes-1,g,g);++fluxReads;
            auto publication=flux.validate();
            for(int face=g-1;face<g+nodes;++face) {
                double ns=(*publication->flux)(mesh.getIdx(face,g,g),publication->massComponent);
                maxFluxDifference=std::max(maxFluxDifference,std::abs(flux.at(0,face,g,g)-ns));++comparedFaces;
            }
            // zeroGradient concentration exterior is the closest physical value.
            auto readC=realized.stageReader("rhoC",context.stageCount),readRho=realized.stageReader("rho",context.stageCount);
            double budget=right*readC(mesh.getIdx(g+nodes-1,g,g),0)/readRho(mesh.getIdx(g+nodes-1,g,g),0)-left*readC(mesh.getIdx(g,g,g),0)/readRho(mesh.getIdx(g,g,g),0);
            net+=(context.stageCount==1?1:system.numericalSystem.time.recipe.finalWeights()[stage]/system.numericalSystem.time.recipe.finalDivisor())*budget;
            log<<step<<','<<stage<<','<<t<<','<<dt<<','<<budget<<'\n';
            c.rhs(context);realized.rhsReady();f.advance(context);c.advance(context);realized.finishStage();
            rejects([&]{flow->stageMassFlux({state.step,stage+1,t});});
            rejects([&]{flux.at(0,g,g,g);});
        }
        f.validatePublish();c.validatePublish();f.publish();c.publish();state.time+=dt;++state.step;accumulated+=dt*net;
        maxBudget=std::max(maxBudget,std::abs(integral()-initial+accumulated));
        for(int i=g;i<g+nodes;++i)if(constant)maxC=std::max(maxC,std::abs(scalar(i,g,g)/mesh(i,g,g,0)-.7));
    }
    require(maxBudget<5e-13,"Shared scalar conservative open-boundary budget failed");require(!constant || maxC<1e-13,"Constant C nonuniform rho compatibility failed");
    require(maxFluxDifference==0,"Borrowed scalar mass flux differs from actual NS face storage");
    std::cout<<"shared "<<FDM::toString(id)<<" constant="<<constant<<" nodes="<<nodes<<" end="<<state.time<<" steps="<<state.step<<" fluxReads="<<fluxReads<<" comparedFaces="<<comparedFaces<<" fluxDifference="<<maxFluxDifference<<" budgetError="<<maxBudget<<" concentrationError="<<maxC<<'\n';
}
void variableDiffusion() {
    const double h=.025;System::FormulaValues v;v.spacing=[=](int axis){return axis==0?h:0;};
    auto rho=[=](int cell,int){return 1+.2*cell*h;};
    auto q=[=](int cell,int){double x=cell*h;return (1+.2*x)*x*x;};
    auto D=[=](int cell,int){double x=cell*h;return (.1+.03*x)/(1+.2*x);};
    v.boundReads={q,rho,[](int,int){return 0.;},D};
    v.boundNeighbor=[&](std::size_t slot,int cell,int axis,int offset,int component){return v.boundReads.at(slot)(cell+offset,component);};
    auto kernel=Discretization::conservativeCentral2Kernel(0,3,false);double error=0;
    for(int cell=1;cell<39;++cell){double x=cell*h;error=std::max(error,std::abs(kernel(v,cell,0)-(2*.1+4*.03*x)));}
    require(error<1e-12,"Variable rho-weighted diffusivity manufactured polynomial failed");
    v.boundReads[1]=[](int,int){return 0.;};rejects([&]{kernel(v,4,0);});v.boundReads[1]=rho;
    v.boundReads[3]=[](int,int){return -1.;};rejects([&]{kernel(v,4,0);});
    std::cout<<"variable Gamma Central2 polynomial error="<<error<<'\n';
}
void multiDimension(int dimensions) {
    const int n=dimensions==3?9:17,g=3,ny=n,nz=dimensions==3?n:1;const double h=1./(n-1),dt=.00025,end=.001;
    FDM::SolverConfig config;config.numerics.maxDeltaT=dt;config.numerics.cfl=2;config.numerics.timeRecipe=FDM::builtInTimeRecipe(FDM::TimeRecipeId::ClassicalRK4);
    auto r=sharedRequest();auto system=System::build(config,r);System::validate(system);
    Field mesh;mesh.setup(n,ny,nz,g,5);auto model=Physics::FluidStateModel::makeSingleFluidPerfectGas(1.4,287.05,0,.72);mesh.setStateModel(model);ScalarField scalar;scalar.setupLike(mesh,"rhoC");
    for(int cell=0;cell<mesh.TotalSize();++cell){int i,j,k;mesh.getIJK(cell,i,j,k);double x=(i-g)*h,y=(j-g)*h,z=(k-g)*h,rho=1+.05*std::sin(2*std::acos(-1.)*x)*std::cos(2*std::acos(-1.)*y);
        mesh.X(i,j,k)=x;mesh.Y(i,j,k)=y;mesh.Z(i,j,k)=z;mesh.Jac(i,j,k)=1/std::pow(h,dimensions);mesh.XiX(i,j,k)=1/h;mesh.EtY(i,j,k)=1/h;mesh.ZeZ(i,j,k)=nz>1?1/h:1;
        mesh(i,j,k,0)=rho;mesh(i,j,k,1)=rho*.1;mesh(i,j,k,2)=rho*.05;mesh(i,j,k,3)=nz>1?rho*.03:0;mesh(i,j,k,4)=2.5+.5*rho*(.01+.0025+(nz>1?.0009:0));scalar(i,j,k)=.7*rho;}
    State::StateBundle state;state.patches={&mesh};state.stateModel=model;state.registerConservativeState();state.distributed.add(State::scalarView("rhoC",0,mesh,scalar));
    auto data=std::make_shared<FDM::ScalarEquationData>();data->target="rhoC";for(auto& bc:data->boundary)bc.kind=FDM::ScalarBoundaryKind::ZeroGradient;
    Execution::LocalRuntime runtime;HomogeneousExterior exterior;FDM::SolverServices services;services.executionRuntime=&runtime;services.boundaryPipeline=&exterior;services.scalarInstances={{"advance.rhoC",data}};
    Time::RunControl control;control.endTime=end;control.writeByStep=true;control.writeIntervalSteps=1000;Application::Execution::executeEquations(config,system,state,services,control);
    double error=0;for(int k=g;k<g+nz;++k)for(int j=g;j<g+ny;++j)for(int i=g;i<g+n;++i)error=std::max(error,std::abs(scalar(i,j,k)/mesh(i,j,k,0)-.7));
    require(error<1e-13,"Multidimensional shared constant concentration compatibility failed");std::cout<<"shared "<<dimensions<<"D nodes="<<n<<"x"<<ny<<"x"<<nz<<" end="<<state.time<<" steps="<<state.step<<" concentration error="<<error<<'\n';
}


double manufactured(int n,double dt,double end,FDM::TimeRecipeId recipe,bool temporal) {
    const int g=3;const double h=1./(n-1),rho=1.2,u=.3,D=temporal?0:.01,pi=std::acos(-1.);
    FDM::SolverConfig config;config.numerics.maxDeltaT=dt;config.numerics.cfl=2;config.numerics.timeRecipe=FDM::builtInTimeRecipe(recipe);
    auto r=sharedRequest();using E=System::FormulaExpr;auto concentration=E::divide(SF::Test::symbol("rhoC"),SF::Test::symbol("rho"));
    auto rhs=E::op("diffusion",{E::multiply(SF::Test::symbol("rho"),E::constantValue(D,Math::signature(0,{0,2,-1,0,0,0,0}))),concentration});
    if(temporal)rhs=E::add(rhs,E::op("source",{SF::Test::symbol("Sc")}));r.userContributions[0].registeredEquations[0].rhs=rhs;
    auto system=System::build(config,r);System::validate(system);
    Field mesh;mesh.setup(n,1,1,g,5);auto model=Physics::FluidStateModel::makeSingleFluidPerfectGas(1.4,287.05,0,.72);mesh.setStateModel(model);ScalarField scalar;scalar.setupLike(mesh,"rhoC");
    const auto exact=[=](double x,double t){return temporal?std::exp(t):1+.1*std::sin(2*pi*(x-u*t))*std::exp(-4*pi*pi*D*t);};
    for(int cell=0;cell<mesh.TotalSize();++cell){int i,j,k;mesh.getIJK(cell,i,j,k);double x=(i-g)*h;mesh.X(i,j,k)=x;mesh.Jac(i,j,k)=1/h;mesh.XiX(i,j,k)=1/h;mesh.EtY(i,j,k)=1;mesh.ZeZ(i,j,k)=1;mesh(i,j,k,0)=rho;mesh(i,j,k,1)=rho*u;mesh(i,j,k,4)=2.5+.5*rho*u*u;scalar(i,j,k)=rho*exact(x,0);}
    State::StateBundle state;state.patches={&mesh};state.stateModel=model;state.registerConservativeState();state.distributed.add(State::scalarView("rhoC",0,mesh,scalar));
    auto data=std::make_shared<FDM::ScalarEquationData>();data->target="rhoC";for(auto& bc:data->boundary)bc.kind=FDM::ScalarBoundaryKind::ZeroGradient;
    for(int side=0;side<2;++side){data->boundary[side].kind=FDM::ScalarBoundaryKind::FixedValue;data->boundary[side].value=[=](double x,double,double,double t){return rho*exact(x,t);};}
    if(temporal)data->sources["Sc"]=[=](double,double,double,double t){return rho*std::exp(t);};
    Execution::LocalRuntime runtime;HomogeneousExterior exterior;FDM::SolverServices services;services.executionRuntime=&runtime;services.boundaryPipeline=&exterior;services.scalarInstances={{"advance.rhoC",data}};
    Time::RunControl control;control.endTime=end;control.writeByStep=true;control.writeIntervalSteps=1000;Application::Execution::executeEquations(config,system,state,services,control);
    double error=0;for(int i=g+1;i<g+n-1;++i)error=std::max(error,std::abs(scalar(i,g,g)/rho-exact((i-g)*h,end)));
    std::cout<<"manufactured shared temporal="<<temporal<<" recipe="<<FDM::toString(recipe)<<" n="<<n<<" dt="<<dt<<" end="<<state.time<<" steps="<<state.step<<" Linf="<<error<<'\n';return error;
}
void convergence() {
    for(auto recipe:{FDM::TimeRecipeId::ForwardEuler,FDM::TimeRecipeId::ClassicalRK4}) {
        double a=manufactured(17,.025,.1,recipe,true),b=manufactured(17,.0125,.1,recipe,true);
        double order=std::log(a/b)/std::log(2.);require(order>(recipe==FDM::TimeRecipeId::ForwardEuler?.9:3.8),"Shared scalar temporal manufactured convergence failed");std::cout<<"shared time order="<<order<<'\n';
    }
    double a=manufactured(33,.0001,.02,FDM::TimeRecipeId::ClassicalRK4,false),b=manufactured(65,.0001,.02,FDM::TimeRecipeId::ClassicalRK4,false),c=manufactured(129,.0001,.02,FDM::TimeRecipeId::ClassicalRK4,false);
    require(a>b && b>c,"Shared scalar spatial refinement did not reduce error");double order=std::log(b/c)/std::log(2.);require(order>.8 && order<1.2,"Shared Upwind1 recipe did not reproduce first-order space convergence");std::cout<<"shared space order="<<order<<'\n';
}

void stageVariableD() {
    const int n=33,g=3;const double rho=1.2,h=1./(n-1),dt=.0005,end=.005;
    FDM::SolverConfig config;config.numerics.maxDeltaT=dt;config.numerics.cfl=2;config.numerics.timeRecipe=FDM::builtInTimeRecipe(FDM::TimeRecipeId::ClassicalRK4);
    auto r=sharedRequest();using E=System::FormulaExpr;r.userContributions[0].registeredEquations[0].rhs=E::op("diffusion",{E::multiply(SF::Test::symbol("rho"),SF::Test::symbol("D")),E::divide(SF::Test::symbol("rhoC"),SF::Test::symbol("rho"))});
    auto coefficient=request("D",.01,false,Math::signature(0,{0,2,-1,0,0,0,0})).userContributions[0];coefficient.records[0].id="user.variableD";coefficient.execution.pop_back();coefficient.execution[0].order=200;r.userContributions.push_back(coefficient);r.composition.solutionVariables.push_back("D");
    auto system=System::build(config,r);System::validate(system);Field mesh;mesh.setup(n,1,1,g,5);auto model=Physics::FluidStateModel::makeSingleFluidPerfectGas(1.4,287.05,0,.72);mesh.setStateModel(model);ScalarField scalar,D;scalar.setupLike(mesh,"rhoC");D.setupLike(mesh,"D");
    for(int cell=0;cell<mesh.TotalSize();++cell){int i,j,k;mesh.getIJK(cell,i,j,k);double x=(i-g)*h;mesh.X(i,j,k)=x;mesh.Jac(i,j,k)=1/h;mesh.XiX(i,j,k)=1/h;mesh.EtY(i,j,k)=1;mesh.ZeZ(i,j,k)=1;mesh(i,j,k,0)=rho;mesh(i,j,k,4)=2.5;scalar(i,j,k)=rho*(1+x);D(i,j,k)=.02+.01*x;}
    State::StateBundle state;state.patches={&mesh};state.stateModel=model;state.registerConservativeState();state.distributed.add(State::scalarView("rhoC",0,mesh,scalar));state.distributed.add(State::scalarView("D",0,mesh,D));
    auto data=std::make_shared<FDM::ScalarEquationData>(),coefs=std::make_shared<FDM::ScalarEquationData>();data->target="rhoC";coefs->target="D";
    for(auto* item:{data.get(),coefs.get()})for(auto& bc:item->boundary)bc.kind=FDM::ScalarBoundaryKind::ZeroGradient;
    for(int side=0;side<2;++side){data->boundary[side].kind=FDM::ScalarBoundaryKind::FixedValue;data->boundary[side].value=[=](double x,double,double,double t){return rho*(1+x+.01*t);};coefs->boundary[side].kind=FDM::ScalarBoundaryKind::FixedValue;coefs->boundary[side].value=[](double x,double,double,double){return .02+.01*x;};}
    Execution::LocalRuntime runtime;HomogeneousExterior exterior;FDM::SolverServices services;services.executionRuntime=&runtime;services.boundaryPipeline=&exterior;services.scalarInstances={{"advance.rhoC",data},{"advance.D",coefs}};
    Time::RunControl control;control.endTime=end;control.writeByStep=true;control.writeIntervalSteps=1000;Application::Execution::executeEquations(config,system,state,services,control);
    double error=0;for(int i=g;i<g+n;++i)error=std::max(error,std::abs(scalar(i,g,g)/rho-(1+(i-g)*h+.01*end)));
    require(error<1e-13,"Variable diffusivity Stage STATE manufactured execution failed");std::cout<<"shared variable Stage D target=rhoC D=.02+.01*x end="<<state.time<<" steps="<<state.step<<" error="<<error<<'\n';
}
void publication() {
    Field f;f.setup(5,1,1,1,5);f.setStateModel(Physics::FluidStateModel::makeSingleFluidPerfectGas(1.4,287.05,0,.72));PatchWorkspace w;w.ensureFor(f);
    rejects([&]{StageMassFluxView(w.massFlux,f,{2,1,.3},"WENO7");});
    // Lifecycle negatives use synthetic values only; actual NS values are tested in shared().
    w.convectiveFlux(f.TotalSize()+f.getIdx(1,1,1),0)=3.25;rejects([&]{w.massFlux->publish(f,w.convectiveFlux,{2,1,.3},"WENO7");});w.massFlux->finalizeCanonical();w.massFlux->publish(f,w.convectiveFlux,{2,1,.3},"WENO7");StageMassFluxView view(w.massFlux,f,{2,1,.3},"WENO7");require(view.at(1,1,1,1)==3.25,"Face orientation/index failed");
    rejects([&]{StageMassFluxView(w.massFlux,f,{2,1,.3},"TENO5");});rejects([&]{view.at(3,1,1,1);});rejects([&]{w.massFlux->publish(f,w.convectiveFlux,{2,1,.3},"WENO7");});w.clear();rejects([&]{view.at(1,1,1,1);});
    w.massFlux->finalizeCanonical();w.massFlux->publish(f,w.convectiveFlux,{3,0,.4},"WENO7");StageMassFluxView moving(w.massFlux,f,{3,0,.4},"WENO7");PatchWorkspace moved(std::move(w));rejects([&]{moving.at(0,1,1,1);});
    Field wrong;wrong.setup(5,1,1,1,4);rejects([&]{wrong.setStateModel(f.stateModel());});moved.massFlux->finalizeCanonical();rejects([&]{moved.massFlux->publish(wrong,moved.convectiveFlux,{3,0,.4},"WENO7");});
}
void missingStageStorage() {
    FDM::SolverConfig config;auto system=System::build(config,sharedRequest());
    Field mesh;mesh.setup(17,1,1,3,5);auto model=Physics::FluidStateModel::makeSingleFluidPerfectGas(1.4,287.05,0,.72);mesh.setStateModel(model);
    for(int cell=0;cell<mesh.TotalSize();++cell){int i,j,k;mesh.getIJK(cell,i,j,k);mesh.X(i,j,k)=(i-3)/16.;}
    ScalarField scalar;scalar.setupLike(mesh,"rhoC");State::StateBundle state;state.patches={&mesh};state.stateModel=model;state.registerConservativeState();state.distributed.add(State::scalarView("rhoC",0,mesh,scalar));
    auto views=system.solvePlan.compiledProgram.stateViews;
    views.erase(std::remove_if(views.begin(),views.end(),[](const auto& view){return view.symbol=="rhoC" && view.kind==System::StateViewKind::Stage;}),views.end());
    auto realized=System::realizeState(system.executableSystem.state,system.runtime,state,views);
    auto data=std::make_shared<FDM::ScalarEquationData>();data->target="rhoC";for(auto& bc:data->boundary)bc.kind=FDM::ScalarBoundaryKind::ZeroGradient;
    FDM::SolverServices services;services.scalarInstances={{"advance.rhoC",data}};double limit=.001;
    // Binding must fail before any numerical flux is requested.
    bool rejected=false,fluxRequested=false;
    try {scalarParticipants(system.numericalSystem,system.solvePlan,state,services,realized,limit,[&](auto)->StageMassFluxView{fluxRequested=true;throw std::runtime_error("Unexpected flux evaluation during missing-STATE binding");});}
    catch(const std::exception& error){std::string message=error.what();std::cout<<"rejected: "<<message<<'\n';rejected=message.find("rhoC")!=std::string::npos && message.find("Provider cannot bind unrequested/owned STATE view")!=std::string::npos;}
    require(rejected && !fluxRequested,"Missing rhoC Stage storage did not fail at STATE binding");
}
}
int main() {
    std::cout.precision(17);publication();missingStageStorage();variableDiffusion();stageVariableD();convergence();multiDimension(2);multiDimension(3);FDM::SolverConfig config;config.numerics.maxDeltaT=.001;
    rejects([&]{System::build(config,sharedRequest(true,true));});
    auto mpi=sharedRequest();mpi.parallel=true;rejects([&]{System::build(config,mpi);});
    auto duplicate=sharedRequest();duplicate.userContributions[0].states.push_back(duplicate.userContributions[0].states[0]);rejects([&]{System::build(config,duplicate);});
    auto duplicateScalar=sharedRequest();duplicateScalar.userContributions[0].execution.insert(duplicateScalar.userContributions[0].execution.begin(),duplicateScalar.userContributions[0].execution[0]);rejects([&]{System::build(config,duplicateScalar);});
    auto duplicateFlux=sharedRequest();System::ExecutionScope duplicateMass;
    duplicateMass.kind=System::ExecutionKind::EquationCall;duplicateMass.id="duplicate.mass";
    duplicateMass.step={"continuity",System::Target{"rho"},"duplicate.mass"};duplicateMass.order=80;
    duplicateFlux.userContributions[0].execution.push_back(duplicateMass);rejects([&]{System::build(config,duplicateFlux);});
    auto unsupportedD=sharedRequest();using E=System::FormulaExpr;
    unsupportedD.userContributions[0].registeredEquations[0].rhs=E::op("diffusion",{E::multiply(SF::Test::symbol("rho"),E::add(E::constantValue(.01,Math::signature(0,{0,2,-1,0,0,0,0})),E::constantValue(.02,Math::signature(0,{0,2,-1,0,0,0,0})))),E::divide(SF::Test::symbol("rhoC"),SF::Test::symbol("rho"))});
    rejects([&]{System::build(config,unsupportedD);});
    auto noFlow=sharedRequest();noFlow.includeDefaultFluidPreset=false;noFlow.singleFluidPreset.reset();rejects([&]{System::build(config,noFlow);});
    for(auto id:{FDM::TimeRecipeId::ForwardEuler,FDM::TimeRecipeId::ClassicalRK4}){shared(id,true);shared(id,false);}
    std::cout<<"shared canonical transport contracts passed\n";
}
