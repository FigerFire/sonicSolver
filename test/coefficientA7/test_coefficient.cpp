#include "solver/system/SF_systemBuilder.h"
#include "solver/system/SF_systemValidator.h"
#include "solver/system/SF_systemPrinter.h"
#include "solver/system/SF_scalarMethod.h"
#include "solver/discretization/diffusion/SF_formulaCentral2.h"
#include "app/application/adapters/SF_scalarResources.h"
#include "core/field/SF_scalarField.h"
#include "core/state/SF_state.h"
#include "app/application/execution/SF_flowLoop.h"
#include "app/application/output/SF_fields.h"
#include "infrastructure/io/SF_caseConfig.h"
#include "infrastructure/io/SF_resultWriter.h"
#include <filesystem>
#include <fstream>
#include <cmath>
#include <iostream>
using namespace SF;using namespace SF::System;
namespace {
using X=FormulaExpr;
void require(bool v,const std::string& why){if(!v)throw std::runtime_error(why);}
template<class F>void rejects(F f,const std::string& expected){try{f();}catch(const std::exception& e){require(std::string(e.what()).find(expected)!=std::string::npos,e.what());std::cout<<"rejected "<<e.what()<<'\n';return;}throw std::runtime_error("Missing rejection: "+expected);}
const auto unit=Math::signature(0,{}),diffDim=Math::signature(0,{0,2,-1,0,0,0,0}),rateDim=Math::signature(0,{0,0,-1,0,0,0,0});
StateSymbol symbol(const std::string& id,Math::Signature signature){StateSymbol s;s.id=s.name=s.storageKey=id;s.signature=signature;s.storageBinding=StorageBinding::NamedDistributed;s.restartEligible=false;return s;}
BuildRequest make(bool dynamic,bool produced,bool reverse=false,bool stage=true){
 BuildRequest request;request.includeDefaultFluidPreset=false;request.composition.stateDeclared=true;request.composition.solutionVariables=dynamic?std::vector<std::string>{"a","b"}:std::vector<std::string>{"a"};SystemContribution c;c.recordContribution("test.coefficient","A7");
 auto a=symbol("a",unit),b=symbol("b",diffDim),D=symbol("D",diffDim),S=symbol("S",rateDim);c.addState(a);
 if(dynamic)c.addState(b);
 if(produced){D.role=StateRole::Derived;D.evaluation=StateEvaluation::Lazy;D.initializationRequired=D.boundaryRequired=D.outputEligible=false;D.dependencies={"b"};auto p=std::make_shared<StateProducer>();p->id="test.twice";p->expression=X::constantValue(2)*b.expression();p->versions={StateViewKind::Physical,StateViewKind::Stage};p->bind=[](Field&,const auto& inputs,State::DistributedFieldView& out)->std::function<void()>{require(inputs.size()==1,"producer binding missing");out.read=[input=inputs[0]](int cell,int){return 2*input.read(cell,0);};return {};};D.producer=p;c.addState(D);}
 else if(!dynamic){D.role=S.role=StateRole::Algebraic;c.addState(D);c.addState(S);}
 auto rhs=X::op(Math::Op::Diffusion,{dynamic&&!produced?b.expression():D.expression(),a.expression()});if(!dynamic)rhs=rhs+S.expression();
 c.addEquation({"a.balance",X::op(Math::Op::Ddt,{a.expression()}),rhs});
 const auto addCall=[&](const std::string& id){ExecutionScope call;call.kind=ExecutionKind::EquationCall;call.id="advance."+id;call.step={id+".balance",{id},call.id};c.addExecution(call);};addCall("a");
 NumericalBinding which{"a.balance","ExplicitCentral2Arithmetic"};if(!dynamic)c.bindInputs({"a.balance",{{"D",StateVersion::Current},{"S",StateVersion::Current}}});else if(produced)c.bindInputs({"a.balance",{{"D",stage?StateVersion::Stage:StateVersion::Current}}});else if(!stage)c.bindInputs({"a.balance",{{"b",StateVersion::Current}}});c.bindNumerics(which);
 if(dynamic){c.addEquation({"b.balance",X::op(Math::Op::Ddt,{b.expression()}),X::constantValue(-.7,rateDim)*b.expression()});addCall("b");c.bindNumerics({"b.balance","ScalarDiffusionCentral2"});}
 if(reverse)std::reverse(c.execution.begin(),c.execution.end());ExecutionScope commit;commit.kind=ExecutionKind::Commit;commit.id="commit";commit.order=1000000;c.addExecution(commit);request.userContributions.push_back(std::move(c));return request;
}
struct Answer{double error=0;std::vector<double> a;};
Answer run(const std::string& name,int n,double dt,double end,bool dynamic,bool produced,bool reverse=false,bool sine=false,bool discrete=false,bool stage=true,bool borrowed=false){
 FDM::SolverConfig config;config.numerics.maxDeltaT=dt;config.numerics.timeRecipe=FDM::builtInTimeRecipe(FDM::TimeRecipeId::ClassicalRK4);auto request=make(dynamic,produced,reverse,stage);if(borrowed)request.userContributions[0].states[1].storageBinding=StorageBinding::ProviderDistributed;auto system=System::build(config,request);System::validate(system);
 Field mesh;mesh.setupGeometry(n,1,1,1);double h=1./(n-1),pi=std::acos(-1.);for(int k=0;k<mesh.MZ();++k)for(int j=0;j<mesh.MY();++j)for(int i=0;i<mesh.MX();++i)mesh.X(i,j,k)=(i-1)*h;
 mesh.getAllSets()["Left"]={mesh.getIdx(1,1,1)};mesh.getAllSets()["Right"]={mesh.getIdx(n,1,1)};
 CaseConfig input;auto description=std::make_shared<Model::Description>();input.modelDescription=description;
 for(const auto& state:system.executableSystem.state.symbols())if(!state.producer){Model::FieldDescriptor f;f.name=f.symbol=state.id;f.initial=0.;if(state.id=="a")for(auto side:{"Left","Right"})f.boundaries[side]={{"type","fixedValue"},{"value",0.}};description->fields.push_back(f);}
 State::StateBundle bundle;bundle.patches={&mesh};ScalarField given;if(borrowed){given.setupLike(mesh,"D",0.);bundle.distributed.add(State::scalarComponentsView("D",0,mesh,{&given},mesh.NG(),State::HaloSyncStage::None));}auto resources=Application::Execution::Detail::bindScalarResources(mesh,bundle,system,input);require(resources.arrays.size()==(borrowed?2:dynamic?2:3),"Coefficient allocation created a second owner");
 for(const auto& state:system.executableSystem.state.symbols())if(!state.producer){auto* field=bundle.distributed.select(state.id,State::HaloSyncStage::None).at(0);for(int i=0;i<n;++i){double x=i*h,val=0;if(state.id=="a")val=(dynamic||sine)?std::sin(pi*x):x*(1-x);if(state.id=="b")val=.1;if(state.id=="D")val=.1+.04*x;if(state.id=="S")val=sine?(.1+.04*x)*pi*pi*std::sin(pi*x)-.04*pi*std::cos(pi*x):.16+.16*x;field->write(mesh.getIdx(i+1,1,1),0,val);}}
 std::filesystem::create_directories("cases/"+name+"/result");std::ofstream("cases/"+name+"/explain.txt")<<describe(system);ResultWriterConfig output;output.caseDir="cases/"+name;output.outputDir="result";output.jobName=name;ResultWriter writer(output);auto fields=Application::Output::registeredStateFields(system.executableSystem.state,bundle,mesh);writer.save(mesh,0,fields);
 FDM::SolverServices services;services.scalarInstances=resources.instances;Time::RunControl control;control.endTime=end;control.writeByStep=true;control.writeIntervalSteps=100000;Application::Execution::executeEquations(config,system,bundle,services,control);writer.save(mesh,bundle.time,fields);
 require(std::abs(bundle.time-end)<1e-12 && mesh.NVar()==0,"Clock/Q regression");Answer result;auto* physical=bundle.distributed.select("a",State::HaloSyncStage::None).at(0);
 const double lambda=discrete?4*std::pow(std::sin(pi*h/2),2)/(h*h):pi*pi;
 double current=1;
 const auto R=[](double z){return 1+z+z*z/2+z*z*z/6+z*z*z*z/24;};
 for(int step=0;step<bundle.step;++step)current*=R(-lambda*(produced?2:1)*.1*std::pow(R(-.7*dt),step)*dt);
 for(int i=0;i<n;++i){double x=i*h,ref=dynamic?std::sin(pi*x)*std::exp(-lambda*(produced?2:1)*.1/.7*(1-std::exp(-.7*end))):sine?std::sin(pi*x):x*(1-x);if(dynamic && !stage)ref=std::sin(pi*x)*current;double value=physical->read(mesh.getIdx(i+1,1,1),0);result.error=std::max(result.error,std::abs(value-ref));result.a.push_back(value);}
 std::ofstream metrics("cases/"+name+"/metrics.json");metrics.precision(17);metrics<<"{\"error\":"<<result.error<<",\"time\":"<<bundle.time<<",\"steps\":"<<bundle.step<<",\"nodes\":"<<n<<",\"dt\":"<<dt<<"}";std::cout<<name<<" error="<<result.error<<" time="<<bundle.time<<" steps="<<bundle.step<<'\n';return result;
}
void faceBudget(){
 const int n=17;const double h=1./(n-1);std::vector<double> D(n),u(n);
 for(int i=0;i<n;++i){double x=i*h;D[i]=.1+.03*std::cos(x);u[i]=x*(1-x)*(1+.2*x);}
 FormulaValues v;v.spacing=[h](int axis){return axis==0?h:0;};
 v.boundReads={[&](int i,int){return D.at(i);},[&](int i,int){return u.at(i);}};
 v.boundNeighbor=[&](std::size_t slot,int i,int,int offset,int){return slot==0?D.at(i+offset):u.at(i+offset);};
 v.coefficientNeighbor=[&](std::size_t slot,int cell,int axis,int offset){return v.boundNeighbor(slot,cell,axis,offset,0);};
 auto kernel=Discretization::central2ArithmeticKernel(0,1,true);double sum=0,power=0,gradient=0;
 for(int i=1;i<n-1;++i){double rhs=kernel(v,i,0);sum+=h*rhs;power+=h*u[i]*rhs;}
 for(int i=0;i<n-1;++i)gradient+=.5*(D[i]+D[i+1])*std::pow(u[i+1]-u[i],2)/h;
 const double boundary=.5*(D[n-2]+D[n-1])*(u[n-1]-u[n-2])/h-.5*(D[0]+D[1])*(u[1]-u[0])/h;
 require(std::abs(sum-boundary)<2e-15,"Internal diffusive faces did not telescope");
 require(std::abs(power+gradient)<2e-15,"Diffusion power/gradient identity");
 std::cout<<"face budget error="<<sum-boundary<<" power error="<<power+gradient<<'\n';
}
void oneWayCross(){
 BuildRequest request;request.includeDefaultFluidPreset=false;request.composition.stateDeclared=true;request.composition.solutionVariables={"a","b"};SystemContribution c;
 auto a=symbol("a",unit),b=symbol("b",unit);c.addState(a);c.addState(b);
 c.addEquation({"a.cross",X::op(Math::Op::Ddt,{a.expression()}),X::op(Math::Op::Diffusion,{X::constantValue(.1,diffDim),b.expression()})});
 c.addEquation({"b.static",X::op(Math::Op::Ddt,{b.expression()}),X::constantValue(0,rateDim)});
 for(auto id:{"a","b"}){ExecutionScope call;call.kind=ExecutionKind::EquationCall;call.id=std::string("advance.")+id;call.step={std::string(id)=="a"?"a.cross":"b.static",{id},call.id};c.addExecution(call);c.bindNumerics({call.step.equation,"ScalarDiffusionCentral2"});}
 ExecutionScope commit;commit.kind=ExecutionKind::Commit;commit.id="commit";commit.order=1000000;c.addExecution(commit);request.userContributions.push_back(std::move(c));
 FDM::SolverConfig config;config.numerics.maxDeltaT=.01;config.numerics.timeRecipe=FDM::builtInTimeRecipe(FDM::TimeRecipeId::ClassicalRK4);auto system=System::build(config,request);System::validate(system);
 Field mesh;const int n=17;const double h=1./(n-1),pi=std::acos(-1.);mesh.setupGeometry(n,1,1,1);
 for(int k=0;k<mesh.MZ();++k)for(int j=0;j<mesh.MY();++j)for(int i=0;i<mesh.MX();++i)mesh.X(i,j,k)=(i-1)*h;
 mesh.getAllSets()["Left"]={mesh.getIdx(1,1,1)};mesh.getAllSets()["Right"]={mesh.getIdx(n,1,1)};
 CaseConfig input;auto description=std::make_shared<Model::Description>();input.modelDescription=description;
 for(auto id:{"a","b"}){Model::FieldDescriptor f;f.name=f.symbol=id;f.initial=0.;for(auto side:{"Left","Right"})f.boundaries[side]={{"type","fixedValue"},{"value",0.}};description->fields.push_back(f);}
 State::StateBundle bundle;bundle.patches={&mesh};auto resources=Application::Execution::Detail::bindScalarResources(mesh,bundle,system,input);
 auto* bdata=bundle.distributed.select("b",State::HaloSyncStage::None).at(0);for(int i=0;i<n;++i)bdata->write(mesh.getIdx(i+1,1,1),0,std::sin(pi*i*h));
 FDM::SolverServices services;services.scalarInstances=resources.instances;Time::RunControl control;control.endTime=.1;control.writeByStep=true;control.writeIntervalSteps=1000;Application::Execution::executeEquations(config,system,bundle,services,control);
 auto* adata=bundle.distributed.select("a",State::HaloSyncStage::None).at(0);double error=0;const double lambda=4*std::pow(std::sin(pi*h/2),2)/(h*h);
 for(int i=0;i<n;++i)error=std::max(error,std::abs(adata->read(mesh.getIdx(i+1,1,1),0)+.1*lambda*bundle.time*std::sin(pi*i*h)));
 require(error<2e-14,"One-way cross diffusion lost target geometry");std::cout<<"oneWayCross error="<<error<<'\n';
 std::filesystem::create_directories("cases/oneWayCross/result");ResultWriterConfig output;output.caseDir="cases/oneWayCross";output.outputDir="result";output.jobName="oneWayCross";ResultWriter writer(output);writer.save(mesh,bundle.time,Application::Output::registeredStateFields(system.executableSystem.state,bundle,mesh));
 std::ofstream("cases/oneWayCross/metrics.json")<<"{\"error\":"<<error<<",\"time\":"<<bundle.time<<",\"steps\":"<<bundle.step<<"}";
}
}
int main(){try{
 std::cout.precision(17);faceBudget();oneWayCross();require(run("borrowedGiven",33,.0001,.05,false,false,false,false,false,true,true).error<2e-13,"Borrowed model coefficient was copied or misbound");require(run("quadratic",33,.0001,.05,false,false).error<2e-13,"div(D grad(a)) omitted coefficient gradient");
 std::vector<double> space;for(int n:{17,33,65})space.push_back(run("space"+std::to_string(n),n,.0001,.05,false,false,false,true).error);require(space[0]/space[1]>3.8 && space[1]/space[2]>3.8,"Variable coefficient spatial order");
 for(bool produced:{false,true}){std::vector<double> error;for(double dt:{.05,.025,.0125})error.push_back(run(std::string(produced?"producer":"stage")+std::to_string(dt),9,dt,.2,true,produced,false,false,true).error);require(error[0]/error[1]>12 && error[1]/error[2]>12,"Coefficient Stage treatment lost RK4 order");auto forward=run(produced?"producerOrder":"stageOrder",9,.0125,.2,true,produced,false,false,true);auto reversed=run(produced?"producerReverse":"stageReverse",9,.0125,.2,true,produced,true,false,true);require(forward.a==reversed.a,"Coefficient read depends on RHS execution order");}
 for(bool produced:{false,true}){
 auto current=run(produced?"producerCurrent":"givenCurrent",9,.05,.2,true,produced,false,false,true,false);
 auto stage=run(produced?"producerStageCompare":"givenStageCompare",9,.05,.2,true,produced,false,false,true,true);
 require(current.error<2e-14,"Current coefficient did not use the physical step-entry value");
 require(std::abs(current.a[4]-stage.a[4])>1e-4,"Current silently refreshed to Stage");}
 FDM::SolverConfig config;auto missing=make(false,false);missing.userContributions[0].callInputs.clear();rejects([&]{System::build(config,missing);},"Missing explicit input version");
 for(auto version:{StateVersion::Frozen,StateVersion::Lagged}){auto wrong=make(false,false);wrong.userContributions[0].callInputs[0].reads[0].version=version;rejects([&]{System::build(config,wrong);},"requires a real snapshot");}
 {auto wrong=make(false,false);wrong.userContributions[0].callInputs[0].reads[0].version=StateVersion::OldTime;rejects([&]{System::build(config,wrong);},"ordinary temporal snapshot");}
 auto unbound=make(false,false);unbound.userContributions[0].callInputs[0].reads[0].version=StateVersion::Stage;rejects([&]{System::build(config,unbound);},"Stage");
 auto unavailable=make(true,true);auto noStage=std::make_shared<StateProducer>(*unavailable.userContributions[0].states.back().producer);noStage->versions={StateViewKind::Physical};unavailable.userContributions[0].states.back().producer=noStage;rejects([&]{System::build(config,unavailable);},"Stage");
 auto conflict=make(true,true);conflict.userContributions[0].callInputs[0].reads.push_back({"b",StateVersion::Current});rejects([&]{System::build(config,conflict);},"dependency version conflict");
 auto producer=make(true,true);auto& p=producer.userContributions[0].states.back();p.evaluation=StateEvaluation::Materialized;rejects([&]{System::build(config,producer);},"explicit materialized Update");
 std::cout<<"A7 coefficient/producer numerical closure passed\n";return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
