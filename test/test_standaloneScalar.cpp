#include "solver/system/SF_systemBuilder.h"
#include "solver/system/SF_systemValidator.h"
#include "solver/system/SF_scalarMethod.h"
#include "solver/system/SF_stateRealizer.h"
#include "core/field/SF_scalarField.h"
#include "app/application/execution/SF_flowLoop.h"
#include "app/application/output/SF_fields.h"
#include "infrastructure/io/SF_resultWriter.h"
#include "infrastructure/execution/SF_localExecutionRuntime.h"
#include "solver/system/SF_systemPrinter.h"
#include "solver/system/SF_providerResolver.h"
#include "core/system/SF_operationIds.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <iostream>
#include <stdexcept>
using namespace SF;
namespace {
void require(bool value,const char* message) {if (!value) throw std::runtime_error(message);}
System::BuildRequest request(const std::string& target,double D,bool source=false) {
    using namespace System;using E=FormulaExpr;
    BuildRequest result;result.includeDefaultFluidPreset=false;
    result.composition.stateDeclared=true;result.composition.solutionVariables={target};
    SystemContribution contribution;
    contribution.recordContribution("user.scalar","standalone scalar equation");
    StateSymbol symbol;symbol.id=symbol.name=target;
    symbol.storageBinding=StorageBinding::NamedDistributed;symbol.storageKey=target;
    symbol.restartEligible=false;
    contribution.addState(symbol);
    auto rhs=E::op("diffusion",{E::constantValue(D),E::symbol(target)});
    if (source) rhs=E::add(std::move(rhs),E::op("source",{E::symbol("Sc")}));
    contribution.addEquation(System::Equation{"transport."+target,E::op("ddt",{E::symbol(target)}),std::move(rhs)});
    ExecutionScope call;call.kind=ExecutionKind::EquationCall;call.id="advance."+target;
    call.step={"transport."+target,{target,TargetKind::Physical},"advance."+target};
    contribution.addExecution(call);
    ExecutionScope commit;commit.kind=ExecutionKind::Commit;commit.id="commit";commit.order=1000000;
    contribution.addExecution(commit);
    contribution.bindNumerics({"transport."+target,"ScalarDiffusionCentral2"});
    result.userContributions.push_back(std::move(contribution));return result;
}
template<class F> void rejects(F&& f) {bool rejected=false;try {f();}catch(const std::exception& e) {rejected=true;std::cout<<"rejected: "<<e.what()<<'\n';}require(rejected,"Expected fail-fast");}
struct Result {double error=0,discreteError=0,time=0,integral=0;int steps=0;std::vector<double> values;};
Result run(const std::string& name,int nodes,double D,double dt,double end,
    FDM::TimeRecipeId method,bool manufactured=false,const std::string& output={},
    bool neumann=false,double amplitude=1,double start=0,bool movingBoundary=false,bool cosine=false) {
    FDM::SolverConfig config;config.numerics.maxDeltaT=dt;
    config.numerics.timeRecipe=FDM::builtInTimeRecipe(method);
    // Scalar execution must not validate unrelated ideal gas/CFL settings.
    config.numerics.cfl=-1;config.numerics.idealGasGamma=0;
    auto system=System::build(config,request(name,D,manufactured));
    Field geometry;geometry.setupGeometry(nodes,1,1,1);
    ScalarField physical;physical.setupLike(geometry,name);
    const double h=1.0/(nodes-1),pi=std::acos(-1.0);
    for (int k=0;k<geometry.MZ();++k)for(int j=0;j<geometry.MY();++j)for(int i=0;i<geometry.MX();++i) {
        const double x=(i-geometry.NG())*h;geometry.X(i,j,k)=x;
        physical(i,j,k)=movingBoundary?std::exp(start)*(1+x)
            :amplitude*(neumann?(cosine?std::cos(pi*x):1):std::sin(pi*x))*std::exp(manufactured?-start:-D*pi*pi*start);
    }
    State::StateBundle state;state.patches={&geometry};state.time=start;state.step=start?17:0;
    state.distributed.add(State::scalarView(name,0,geometry,physical));
    FDM::ScalarEquationData data;data.target=name;
    for (auto& bc:data.boundary) bc.kind=FDM::ScalarBoundaryKind::ZeroGradient;
    if (!neumann) for (int side=0;side<2;++side) {
        data.boundary[side].kind=FDM::ScalarBoundaryKind::FixedValue;
        data.boundary[side].value=movingBoundary?FDM::ScalarKnownValue([](double x,double,double,double t) {return std::exp(t)*(1+x);})
            :FDM::ScalarKnownValue([](double,double,double,double) {return 0;});
    }
    std::vector<double> sourceTimes;
    if (manufactured) data.sources["Sc"]=[&](double x,double,double,double t) {
        sourceTimes.push_back(t);
        return movingBoundary?std::exp(t)*(1+x):(D*pi*pi-1)*amplitude*std::exp(-t)*std::sin(pi*x);
    };
    Execution::LocalRuntime runtime;FDM::SolverServices services;
    services.executionRuntime=&runtime;services.scalarEquation=&data;
    Time::RunControl control;control.startTime=start;control.endTime=end;control.writeIntervalTime=(end-start)/2;
    if (output.empty()) {control.writeByStep=true;control.writeIntervalSteps=100000;}
    ResultWriterConfig writeConfig;writeConfig.caseDir=output;writeConfig.outputDir="result";writeConfig.jobName=name;
    ResultWriter writer(writeConfig);
    const auto fields=Application::Output::registeredStateFields(system.executableSystem.state,state,geometry);
    require(fields.size()==1 && fields.front().name==name,"Physical STATE output not registered");
    if (!output.empty()) {std::filesystem::create_directories(output+"/result");writer.save(geometry,start,fields);
        std::ofstream(output+"/explain.txt")<<System::describe(system);}
    std::vector<double> dts;
    Time::DriverCallbacks callbacks;
    callbacks.saveTime=[&](double t) {if (!output.empty()) writer.save(geometry,t,fields);};
    callbacks.report=[&](int,const Time::AdvanceResult& result) {dts.push_back(result.dt);};
    Application::Execution::executeEquations(config,system,state,services,control,callbacks);
    require(geometry.NVar()==0 && !state.stateModel && !geometry.stateModel(),"Scalar execution acquired Q/EOS");
    require(std::abs(state.time-end)<1e-12,"Scalar did not reach endTime");
    require(state.step==static_cast<int>(dts.size())+(start?17:0),"Scalar physical step clock reset/duplicated");
    Result result;result.time=state.time;result.steps=state.step;
    const double discreteLambda=4*D*std::pow(std::sin(pi*h/2),2)/(h*h);
    for (int i=0;i<nodes;++i) {
        const double x=i*h,value=physical(i+1,1,1);
        const double exact=neumann?(cosine?amplitude*std::cos(pi*x)*std::exp(-D*pi*pi*end):amplitude):movingBoundary?std::exp(end)*(1+x)
            :amplitude*(neumann?1:std::sin(pi*x))*std::exp(manufactured?-end:-D*pi*pi*end);
        const double discrete=amplitude*std::sin(pi*x)*std::exp(-D*pi*pi*start-discreteLambda*(end-start));
        result.error=std::max(result.error,std::abs(value-exact));
        if (!manufactured && !neumann && !movingBoundary)
            result.discreteError=std::max(result.discreteError,std::abs(value-discrete));
        result.integral+=value*h*(i==0 || i+1==nodes?.5:1);
        result.values.push_back(value);
    }
    if (manufactured) {
        require(!sourceTimes.empty(),"Manufactured source never evaluated");
        if (method==FDM::TimeRecipeId::ClassicalRK4) {
            const int interior=nodes-2;
            require(std::abs(sourceTimes.at(interior)- (start+.5*dts.front()))<1e-12,"Source ignored RK4 stage time");
            require(std::abs(sourceTimes.at(3*interior)-(start+dts.front()))<1e-12,"Source stage4 time wrong");
        }
    }
    std::cout<<"scalar "<<name<<" nodes="<<nodes<<" dt="<<dt<<" end="<<end
        <<" error="<<result.error<<" timeError="<<result.discreteError<<" integral="<<result.integral<<'\n';
    if (!output.empty()) {
        std::ofstream metrics(output+"/metrics.txt");metrics.precision(17);
        metrics<<"time "<<state.time<<"\nsteps "<<state.step<<"\nerror "<<result.error<<"\nintegral "<<result.integral<<'\n';
        for (double value:dts) metrics<<"dt "<<value<<'\n';
    }
    return result;
}
void multidimensional(int dimensions) {
    const int n=17;const double D=.03,end=.05,h=1./(n-1),pi=std::acos(-1.);
    const std::string name="scalar"+std::to_string(dimensions)+"D";
    FDM::SolverConfig config;config.numerics.maxDeltaT=.001;
    config.numerics.timeRecipe=FDM::builtInTimeRecipe(FDM::TimeRecipeId::ClassicalRK4);
    auto system=System::build(config,request(name,D));
    Field mesh;mesh.setupGeometry(n,n,dimensions==3?n:1,1);
    ScalarField scalar;scalar.setupLike(mesh,name);
    for (int k=0;k<mesh.MZ();++k)for(int j=0;j<mesh.MY();++j)for(int i=0;i<mesh.MX();++i) {
        const double x=(i-1)*h,y=(j-1)*h,z=dimensions==3?(k-1)*h:0;
        mesh.X(i,j,k)=x;mesh.Y(i,j,k)=y;mesh.Z(i,j,k)=z;
        scalar(i,j,k)=std::sin(pi*x)*std::sin(pi*y)*(dimensions==3?std::sin(pi*z):1);
    }
    State::StateBundle state;state.patches={&mesh};state.distributed.add(State::scalarView(name,0,mesh,scalar));
    FDM::ScalarEquationData data;data.target=name;
    for (auto& bc:data.boundary) {bc.kind=FDM::ScalarBoundaryKind::FixedValue;bc.value=[](double,double,double,double) {return 0;};}
    FDM::SolverServices services;services.scalarEquation=&data;
    Time::RunControl control;control.endTime=end;control.writeByStep=true;control.writeIntervalSteps=1000;
    Application::Execution::executeEquations(config,system,state,services,control);
    double error=0;
    for (int k=1;k<mesh.NZ()+1;++k)for(int j=1;j<n+1;++j)for(int i=1;i<n+1;++i) {
        const double exact=std::sin(pi*(i-1)*h)*std::sin(pi*(j-1)*h)
            *(dimensions==3?std::sin(pi*(k-1)*h):1)*std::exp(-dimensions*D*pi*pi*end);
        error=std::max(error,std::abs(scalar(i,j,k)-exact));
    }
    require(error<.00015,"Multidimensional scalar diffusion failed");
    std::cout<<dimensions<<"D analytic error "<<error<<'\n';
}
}
int main(int argc,char** argv) {
    std::cout.precision(17);
    FDM::SolverConfig config;config.numerics.maxDeltaT=.001;
    for (const auto& name:{"c","tracer","scalarA","scalarB"}) {
        auto system=System::build(config,request(name,.01,true));System::validate(system);
        require(system.rawSystem.state.size()==1,"Scalar composition inherited flow STATE");
        require(!system.executableSystem.registry.contains("momentum"),"Scalar inherited Momentum");
        require(system.runtime.report.status==System::RuntimeStatus::Runnable,"Scalar plan unresolved");
        for (const auto& op:system.runtime.operationBindings)
            require(op.provider==System::ScalarOps::Provider || op.provider==System::TemporalOps::Provider,"Scalar acquired flow provider");
    }
    Field geometry;geometry.setupGeometry(11,1,1,1);
    require(geometry.NVar()==0 && !geometry.stateModel(),"Scalar geometry allocated Q/EOS");
    ScalarField scalar;scalar.setupLike(geometry,"tracer");
    State::StateBundle state;state.patches={&geometry};state.distributed.add(State::scalarView("tracer",0,geometry,scalar));
    auto system=System::build(config,request("tracer",.01));
    auto realized=System::realizeState(system.executableSystem.state,system.runtime,state,system.solvePlan.compiledProgram.stateViews);
    require(realized.at("tracer").fields.size()==1,"Scalar physical storage missing");
    rejects([&] {state.distributed.add(State::scalarView("tracer",0,geometry,scalar));});
    auto absent=request("tracer",.01);absent.userContributions[0].states.clear();rejects([&] {System::build(config,absent);});
    auto missing=request("tracer",.01);missing.userContributions[0].numerics.clear();rejects([&] {System::build(config,missing);});
    auto unknown=request("tracer",.01);unknown.userContributions[0].numerics[0].method="ImplicitUnknown";rejects([&] {System::build(config,unknown);});
    auto mpi=request("tracer",.01);mpi.parallel=true;rejects([&] {System::build(config,mpi);});
    rejects([&] {System::build(config,request("tracer",-1));});
    auto bad=config;bad.numerics.maxDeltaT=0;rejects([&] {System::build(bad,request("tracer",.01));});
    for (double dt:{-1.0,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}) {
        auto invalid=config;invalid.numerics.maxDeltaT=dt;
        rejects([&] {System::build(invalid,request("tracer",.01));});
    }
    auto unsupported=request("tracer",.01);
    unsupported.userContributions[0].registeredEquations[0].rhs=System::FormulaExpr::op("div",{System::FormulaExpr::symbol("tracer")});
    rejects([&] {System::build(config,unsupported);});
    auto vector=request("tracer",.01);vector.userContributions[0].states[0].shape=System::ValueShape::Vector;
    vector.userContributions[0].states[0].components=3;rejects([&] {System::build(config,vector);});
    auto working=request("tracer",.01);working.userContributions[0].execution[0].step.target.kind=System::TargetKind::Working;
    rejects([&] {System::build(config,working);});
    auto duplicate=request("tracer",.01);duplicate.userContributions[0].states.push_back(duplicate.userContributions[0].states[0]);
    rejects([&] {System::build(config,duplicate);});
    auto missingCommit=request("tracer",.01);missingCommit.userContributions[0].execution.pop_back();
    rejects([&] {System::build(config,missingCommit);});
    auto ssp=config;ssp.numerics.timeRecipe=FDM::builtInTimeRecipe(FDM::TimeRecipeId::SSPRK3);
    rejects([&] {System::build(ssp,request("tracer",.01));});
    auto implicit=request("tracer",.01);implicit.userContributions[0].numerics[0].method="LinearEquation";
    rejects([&] {System::build(config,implicit);});
    auto symbolic=request("tracer",.01);
    System::StateSymbol coefficient;coefficient.id=coefficient.name="materialD";
    coefficient.role=System::StateRole::Algebraic;coefficient.constantValue=.01;
    coefficient.runtimeStorageRequired=false;coefficient.initializationRequired=false;coefficient.boundaryRequired=false;
    symbolic.userContributions[0].addState(coefficient);
    symbolic.userContributions[0].registeredEquations[0].rhs.arguments[0]=System::FormulaExpr::symbol("materialD");
    auto symbolicSystem=System::build(config,symbolic);System::validate(symbolicSystem);
    const auto* symbolicContract=std::any_cast<System::CompiledScalarEquation>(&symbolicSystem.solvePlan.compiledProgram.steps[0].providerContract);
    require(symbolicContract && symbolicContract->diffusivity==.01,"Constant STATE diffusivity not compiled");
    auto corrupted=system;
    const auto corrupt=[&](const auto& self,System::SolvePlanNode& node)->void {
        if (node.operation==system.solvePlan.compiledProgram.steps.front().temporalRhs) node.operation=System::OpIds::ExplicitStageExecute;
        for (auto& child:node.children) self(self,child);
    };
    corrupt(corrupt,corrupted.solvePlan.root);
    auto corruptBindings=System::compileOperationBindings(corrupted.executableSystem,corrupted.numericalSystem,corrupted.solvePlan,{},false);
    require(System::reportOperationBindings(corrupted.solvePlan,corruptBindings).status!=System::RuntimeStatus::Runnable,"Flow-only operation falsely Runnable");
    // Missing real storage/boundary/source and wrong geometry must fail during
    // stable binding, before the first numerical stage or physical clock commit.
    const auto invalidBinding=[&](int mode) {
        Field mesh;mesh.setupGeometry(9,1,1,1);
        for (int k=0;k<mesh.MZ();++k)for(int j=0;j<mesh.MY();++j)for(int i=0;i<mesh.MX();++i) mesh.X(i,j,k)=(i-1)/8.;
        ScalarField value;value.setupLike(mesh,"tracer",1);
        State::StateBundle bundle;bundle.patches={&mesh};
        if (mode!=0) bundle.distributed.add(State::scalarView("tracer",0,mesh,value));
        if (mode==3) bundle.patches.push_back(&mesh);
        if (mode==4) mesh.X(4,1,1)=.4;
        FDM::ScalarEquationData data;data.target=mode==5?"other":"tracer";
        for (auto& bc:data.boundary) bc.kind=FDM::ScalarBoundaryKind::ZeroGradient;
        if (mode==1) data.boundary[0].kind=FDM::ScalarBoundaryKind::Unspecified;
        FDM::SolverServices services;services.scalarEquation=&data;
        auto resolved=System::build(config,request("tracer",.01,mode==2));
        Time::RunControl control;control.endTime=.01;control.writeByStep=true;control.writeIntervalSteps=1;
        rejects([&] {Application::Execution::executeEquations(config,resolved,bundle,services,control);});
        require(bundle.time==0 && bundle.step==0,"Rejected binding committed a physical step");
    };
    for (int mode=0;mode<6;++mode) invalidBinding(mode);
    std::cout<<"standalone scalar composition passed\n";
    if (argc>1 && std::string(argv[1])=="--composition") return 0;
    const auto euler=FDM::TimeRecipeId::ForwardEuler,rk4=FDM::TimeRecipeId::ClassicalRK4;
    const std::string root=argc>1?argv[1]:"standalone-scalar-output";
    require(run("c",65,.01,.001,1,euler,false,root+"/euler").error<3e-5,"Euler analytic error");
    require(run("tracer",65,.01,.001,1,rk4,false,root+"/rk4").error<3e-5,"RK4 analytic error");
    const auto a=run("scalarA",33,.03,.001,.1,rk4,true);
    require(a.error<1e-4,"RK4 manufactured source error");
    require(run("sourceEuler",33,.03,.0005,.1,euler,true).error<1e-4,"Euler source manufactured solution failed");
    const auto renamed=run("scalarB",33,.03,.001,.1,rk4,true);
    require(a.values==renamed.values,"Renamed STATE/equation changed numerical result");
    const auto independent=run("independent",25,.08,.0005,.1,rk4,true,{},false,2);
    require(independent.values!=a.values && independent.error<.001,"Second independent system failed");
    require(run("constant",33,.03,.001,.1,rk4,false,{},true,2).error<1e-13,"zeroGradient constant solution changed");
    require(run("cosine",33,.03,.001,.1,rk4,false,{},true,1,0,false,true).error<3e-5,"zeroGradient diffusion boundary failed");
    require(run("movingDirichlet",17,.01,.001,.1,rk4,true,{},false,1,0,true).error<1e-10,"Stage-time boundary/source mismatch");
    require(run("restart",33,.03,.001,.31,rk4,false,{},false,1,.3).error<1e-4,"Restart clock/state failed");
    multidimensional(2);multidimensional(3);
    std::vector<double> space;
    for (int n:{17,33,65}) space.push_back(run("space",n,.1,.0001,.1,rk4).error);
    require(space[0]/space[1]>3.9 && space[1]/space[2]>3.9,"Central2 spatial refinement failed");
    for (auto method:{euler,rk4}) {
        std::vector<double> errors;
        for (double dt:{.02,.01,.005}) errors.push_back(run("time",9,.1,dt,.1,method).discreteError);
        const double expected=method==euler?1.9:14;
        require(errors[0]/errors[1]>expected && errors[1]/errors[2]>expected,"Temporal order against semi-discrete reference failed");
        std::cout<<"temporal ratios "<<errors[0]/errors[1]<<' '<<errors[1]/errors[2]<<'\n';
    }
    std::cout<<"standalone scalar execution passed\n";
    return 0;
}
