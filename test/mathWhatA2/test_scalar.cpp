#include "solver/system/SF_systemBuilder.h"
#include "solver/system/SF_systemValidator.h"
#include "solver/system/SF_systemPrinter.h"
#include "solver/system/SF_scalarMethod.h"
#include "app/application/execution/SF_flowLoop.h"
#include "app/application/output/SF_fields.h"
#include "infrastructure/io/SF_resultWriter.h"
#include "infrastructure/execution/SF_localExecutionRuntime.h"
#include "core/field/SF_scalarField.h"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
using namespace SF;
namespace {
void require(bool condition,const std::string& message){if(!condition)throw std::runtime_error(message);}
struct Spec {std::string id,target="c";FDM::TimeRecipeId recipe=FDM::TimeRecipeId::ClassicalRK4;double dt=.05,end=.5;int nodes=9;bool power=false,extra=false,manufactured=false,leftDiffusion=false;};
struct Result {std::vector<double> values;double error=0;};
Result run(const std::filesystem::path& base,const Spec& spec) {
    using namespace System;using E=FormulaExpr;
    FDM::SolverConfig config;config.numerics.timeRecipe=FDM::builtInTimeRecipe(spec.recipe);config.numerics.maxDeltaT=spec.dt;
    BuildRequest request;request.includeDefaultFluidPreset=false;request.composition.stateDeclared=true;
    SystemContribution contribution;contribution.recordContribution("test.reactionDiffusion","typed explicit scalar composition");
    std::vector<std::string> targets{spec.target};if(spec.extra)targets.push_back("independent");
    for(const auto& target:targets) {
        StateSymbol symbol;symbol.id=symbol.name=symbol.storageKey=target;symbol.storageBinding=StorageBinding::NamedDistributed;symbol.signature=Math::signature(0,{});symbol.restartEligible=false;
        contribution.addState(symbol);request.composition.solutionVariables.push_back(target);const auto c=symbol.expression();
        const auto k=E::constantValue(target==spec.target?.7:1.1,Math::signature(0,{0,0,-1,0,0,0,0}));
        const auto diffusion=E::op(Math::Op::Diffusion,{E::constantValue(.03,Math::signature(0,{0,2,-1,0,0,0,0})),c},"diffusion."+target);
        auto reaction=spec.power?E::negate(k*E::power(c,2)):E::negate(k)*c*c;
        auto rhs=spec.leftDiffusion?reaction:spec.power?diffusion-k*E::power(c,2):diffusion+reaction;
        if(spec.manufactured)rhs=rhs+E::op(Math::Op::Source,{E::symbol("manufactured",Math::signature(0,{0,0,-1,0,0,0,0}))},"source."+target);
        const auto transient=E::op(Math::Op::Ddt,{c},"time."+target);
        contribution.addEquation({"equation."+target,spec.leftDiffusion?transient-diffusion:transient,rhs,{OriginKind::User,"A2 "+spec.id}});
        ExecutionScope call;call.kind=ExecutionKind::EquationCall;call.step={"equation."+target,{target},"advance."+target};contribution.addExecution(call);
        contribution.bindNumerics({"equation."+target,"ScalarDiffusionCentral2"});
    }
    ExecutionScope commit;commit.kind=ExecutionKind::Commit;commit.order=100;commit.id="physical.commit";contribution.addExecution(commit);request.userContributions.push_back(contribution);
    auto system=System::build(config,request);System::validate(system);
    const auto path=base/spec.id;std::filesystem::create_directories(path/"result");
    std::ofstream(path/"explain.txt")<<System::describe(system);
    std::ofstream settings(path/"settings.txt");settings<<std::setprecision(17)<<"equation ddt(c)=diffusion(0.03,c)-0.7*c^2; no flow/EOS/models/IBM/MPI\nSTATE dimensionless scalar NamedDistributed; old/stage/RHS original scalar workspace\nx=[0,1]; y,z inactive; nodes="<<spec.nodes<<" NG=1\nIC="<<(spec.manufactured?"exp(-t)*(1+0.2*cos(pi*x)); manufactured source":"uniform c0=1; independent c0=1.3, k=1.1 if selected")<<"\nBC six zeroGradient; Central2; start=0 step=0 end="<<spec.end<<" maxDeltaT="<<spec.dt<<" recipe="<<FDM::toString(spec.recipe)<<"\npower="<<spec.power<<" extra="<<spec.extra<<" leftDiffusion="<<spec.leftDiffusion<<"\n";
    Field geometry;geometry.setupGeometry(spec.nodes,1,1,1);const double h=1./(spec.nodes-1),pi=std::acos(-1.);
    for(int z=0;z<geometry.MZ();++z)for(int y=0;y<geometry.MY();++y)for(int x=0;x<geometry.MX();++x)geometry.X(x,y,z)=(x-geometry.NG())*h;
    State::StateBundle state;state.patches={&geometry};std::vector<std::unique_ptr<ScalarField>> fields;
    std::vector<std::shared_ptr<FDM::ScalarEquationData>> retained;FDM::SolverServices services;Execution::LocalRuntime runtime;services.executionRuntime=&runtime;
    for(const auto& target:targets) {
        auto field=std::make_unique<ScalarField>();field->setupLike(geometry,target);
        for(int z=0;z<geometry.MZ();++z)for(int y=0;y<geometry.MY();++y)for(int x=0;x<geometry.MX();++x)(*field)(x,y,z)=spec.manufactured?1+.2*std::cos(pi*geometry.X(x,y,z)):target==spec.target?1:1.3;
        state.distributed.add(State::scalarView(target,0,geometry,*field));fields.push_back(std::move(field));
        auto data=std::make_shared<FDM::ScalarEquationData>();data->target=target;for(auto& bc:data->boundary)bc.kind=FDM::ScalarBoundaryKind::ZeroGradient;
        if(spec.manufactured)data->sources["manufactured"]=[pi](double x,double,double,double t){const double c=std::exp(-t)*(1+.2*std::cos(pi*x));return -c+.03*pi*pi*.2*std::cos(pi*x)*std::exp(-t)+.7*c*c;};
        services.scalarInstances.push_back({"advance."+target,data});retained.push_back(std::move(data));
    }
    ResultWriterConfig output;output.caseDir=path.string();output.outputDir="result";output.jobName="scalar";ResultWriter writer(output);
    const auto outputFields=Application::Output::registeredStateFields(system.executableSystem.state,state,geometry);writer.save(geometry,0,outputFields);
    Time::RunControl control;control.endTime=spec.end;control.writeByStep=true;control.writeIntervalSteps=100000;
    std::ofstream times(path/"dt.txt");times<<std::setprecision(17);Time::DriverCallbacks callbacks;callbacks.report=[&](int,const Time::AdvanceResult& r){times<<r.time<<' '<<r.dt<<'\n';};
    Application::Execution::executeEquations(config,system,state,services,control,callbacks);writer.save(geometry,state.time,outputFields);
    require(std::abs(state.time-spec.end)<1e-14,"Physical clock not committed correctly");
    Result result;std::ofstream values(path/"result/state.csv");values<<std::setprecision(17)<<"x,c,exact,error\n";
    for(int i=0;i<spec.nodes;++i) {
        const double x=i*h,value=(*fields.front())(i+1,1,1),exact=spec.manufactured?std::exp(-spec.end)*(1+.2*std::cos(pi*x)):1/(1+.7*spec.end);
        result.values.push_back(value);result.error=std::max(result.error,std::abs(value-exact));values<<x<<','<<value<<','<<exact<<','<<value-exact<<'\n';
    }
    std::ofstream metrics(path/"metrics.txt");metrics<<std::setprecision(17)<<"time "<<state.time<<"\nsteps "<<state.step<<"\nLinf "<<result.error<<'\n';
    std::cout<<spec.id<<" steps="<<state.step<<" final="<<state.time<<" Linf="<<std::setprecision(17)<<result.error<<'\n';return result;
}
}
int main(int argc,char** argv) {
 try {
    const std::filesystem::path base=argc>1?argv[1]:"cases";std::filesystem::create_directories(base);
    for(auto method:{FDM::TimeRecipeId::ForwardEuler,FDM::TimeRecipeId::ClassicalRK4}) {
        std::vector<double> error;const auto label=method==FDM::TimeRecipeId::ForwardEuler?"euler":"rk4";
        for(double dt:{.1,.05,.025}) {Spec s;s.id=std::string(label)+"-"+std::to_string(dt);s.recipe=method;s.dt=dt;error.push_back(run(base,s).error);}
        std::cout<<label<<" error ratios "<<error[0]/error[1]<<' '<<error[1]/error[2]<<'\n';
        require(error[0]/error[1]>(method==FDM::TimeRecipeId::ForwardEuler?1.8:12),"Unexpected first temporal convergence ratio");
        require(error[1]/error[2]>(method==FDM::TimeRecipeId::ForwardEuler?1.8:12),"Unexpected second temporal convergence ratio");
    }
    Spec reference;reference.id="multiply";auto normal=run(base,reference);
    auto changed=reference;changed.id="power";changed.power=true;auto power=run(base,changed);
    changed.id="renamed";changed.target="tracer";auto renamed=run(base,changed);
    changed.id="independent-instance";changed.extra=true;auto extra=run(base,changed);
    changed=reference;changed.id="left-diffusion";changed.leftDiffusion=true;auto left=run(base,changed);
    double difference=0;for(std::size_t i=0;i<normal.values.size();++i) {
        difference=std::max(difference,std::abs(normal.values[i]-power.values[i]));
        require(power.values[i]==renamed.values[i] && renamed.values[i]==extra.values[i],"Names or an independent equation contaminated the original solution");
        require(std::abs(normal.values[i]-left.values[i])<3e-14,"Controlled residual extraction changed the scalar meaning");
    }
    require(difference<3e-14,"Multiply/Pow arithmetic differ beyond floating precision");std::cout<<"multiply/power maxDifference="<<difference<<'\n';
    std::vector<double> spatial;
    for(int n:{17,33}) {Spec s;s.id="nonuniform-"+std::to_string(n);s.nodes=n;s.dt=.0001;s.end=.02;s.power=true;s.manufactured=true;spatial.push_back(run(base,s).error);}
    require(spatial[0]/spatial[1]>3.5 && spatial[1]<1e-5,"Nonuniform diffusion+nonlinear reaction/source failed independent manufactured reference");
    std::cout<<"manufactured spatial ratio "<<spatial[0]/spatial[1]<<"\nA2 nonlinear scalar execution passed\n";
 } catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
