#include "solver/system/SF_systemBuilder.h"
#include "solver/system/SF_systemValidator.h"
#include "solver/system/SF_scalarMethod.h"
#include "app/application/adapters/SF_scalarResources.h"
#include "solver/algorithm/SF_singleFluidStepper.h"
#include "core/field/SF_scalarField.h"
#include "core/state/SF_state.h"
#include "infrastructure/io/SF_caseConfig.h"
#include "infrastructure/mpi/backend/SF_mpiBackend.h"
#include "app/application/execution/SF_flowLoop.h"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <iostream>
using namespace SF;using namespace SF::System;
namespace {
using X=FormulaExpr;
void require(bool value,const std::string& message){if(!value)throw std::runtime_error(message);}
template<class F>void rejects(F function,const std::string& part){try{function();}catch(const std::exception& e){require(std::string(e.what()).find(part)!=std::string::npos,e.what());std::cout<<"rejected "<<e.what()<<'\n';return;}throw std::runtime_error("Missing rejection "+part);}
const auto unit=Math::signature(0,{}),rate=Math::signature(0,{0,0,-1,0,0,0,0}),diff=Math::signature(0,{0,2,-1,0,0,0,0});
StateSymbol symbol(std::string name,Math::Signature dim){StateSymbol s;s.id=s.name=s.storageKey=name;s.signature=dim;s.storageBinding=StorageBinding::NamedDistributed;s.restartEligible=false;return s;}
BuildRequest request(bool badLayout=false,bool badVersion=false,bool nonconverged=false){
 BuildRequest r;r.includeDefaultFluidPreset=false;r.composition.stateDeclared=true;r.composition.solutionVariables={"a","b"};SystemContribution c;c.recordContribution("test.block","A8");
 auto a=symbol("a",unit),b=symbol("b",unit),f=symbol("f",rate),g=symbol("g",rate);f.role=StateRole::Algebraic;f.storageBinding=StorageBinding::ProviderDistributed;if(badLayout)f.componentOffset=1;c.addState(a);c.addState(b);c.addState(f);c.addState(g);
 const auto diffusion=[&](auto s,double D){return X::op(Math::Op::Diffusion,{X::constantValue(D,diff),s.expression()});};
 c.addEquation({"G",-diffusion(a,.1)+X::constantValue(2,rate)*a.expression()-X::constantValue(.4,rate)*b.expression(),f.expression()});
 c.addEquation({"B",-diffusion(b,.2)+X::constantValue(.3,rate)*a.expression()+X::constantValue(1.5,rate)*b.expression(),g.expression()});
 ExecutionScope block;block.kind=ExecutionKind::BlockSolve;block.id="whole";block.equations={"G","B"};block.unknowns={"a","b"};c.addExecution(block);
 NumericalBinding which{"*","LinearBlockCentral2Hypre"};which.occurrence="whole";which.parameters={{"relativeTolerance",nonconverged?1e-15:1e-12},{"absoluteTolerance",nonconverged?0:1e-13},{"maxIterations",nonconverged?1:500}};if(badVersion)c.bindInputs({"whole",{{"f",StateVersion::Lagged}},{},ReferenceKind::Step});c.bindNumerics(which);r.userContributions.push_back(c);return r;
}
void run(bool expire,bool fail){
 FDM::SolverConfig config;config.numerics.timeRecipe=FDM::builtInTimeRecipe(FDM::TimeRecipeId::ClassicalRK4);auto system=System::build(config,request(false,false,fail));System::validate(system);
 const int n=17;Field mesh;mesh.setupGeometry(n,1,1,1);for(int k=0;k<mesh.MZ();++k)for(int j=0;j<mesh.MY();++j)for(int i=0;i<mesh.MX();++i)mesh.X(i,j,k)=(i-1.)/(n-1);
 mesh.getAllSets()["Left"]={mesh.getIdx(1,1,1)};mesh.getAllSets()["Right"]={mesh.getIdx(n,1,1)};
 State::StateBundle state;state.patches={&mesh};state.time=3.;state.step=7;state.dt=.125;
 auto original=std::make_shared<ScalarField>();original->setupLike(mesh,"f",0);std::weak_ptr<ScalarField> weak=original;
 for(int i=0;i<n;++i){double x=i/double(n-1),a=1+x*(1-x),b=2+.5*x*(1-x);original->values()[mesh.getIdx(i+1,1,1)]=.2+2*a-.4*b;}
 auto borrowed=State::scalarComponentsView("f",0,mesh,{original.get()},1,State::HaloSyncStage::None);borrowed.requireReadable=[weak]{if(weak.expired())throw std::runtime_error("External model STATE owner expired: f");};state.distributed.add(borrowed);
 CaseConfig input;auto description=std::make_shared<Model::Description>();input.modelDescription=description;
 for(const auto& name:{"a","b","g"}){Model::FieldDescriptor f;f.name=f.symbol=name;f.initial=name==std::string("a")?.3:name==std::string("b")?.4:0.;if(name!=std::string("g"))for(const auto& side:{"Left","Right"})f.boundaries[side]={{"type","fixedValue"},{"value",name==std::string("a")?1.:2.}};description->fields.push_back(f);}
 auto resources=Application::Execution::Detail::bindScalarResources(mesh,state,system,input);require(resources.arrays.size()==3,"External model resource copied to second physical array");
 auto g=state.distributed.select("g",State::HaloSyncStage::None).at(0);for(int i=0;i<n;++i){double x=i/double(n-1),a=1+x*(1-x),b=2+.5*x*(1-x);g->write(mesh.getIdx(i+1,1,1),0,.2+.3*a+1.5*b);}
 FDM::SolverServices services;services.scalarInstances=resources.instances;SolverAlgorithm::SingleFluidStepper host(config,system.executableSystem,system.numericalSystem,system.solvePlan,system.runtime);host.bindServices(services);host.bindSolvePlan(system.solvePlan);FDM::SolverState active;active.bundle=&state;host.prepare(active);
 if(expire)original.reset();
 if(expire || fail){rejects([&]{host.advance(active);},expire?"owner expired":"did not converge");for(const auto& target:{"a","b"}){auto field=state.distributed.select(target,State::HaloSyncStage::None).at(0);for(int i=0;i<n;++i)require(field->read(mesh.getIdx(i+1,1,1),0)==(target==std::string("a")?.3:.4),"Failed block wrote a partial physical solution");}}
 else {host.advance(active);double error=0;for(const auto& target:{"a","b"}){auto field=state.distributed.select(target,State::HaloSyncStage::None).at(0);for(int i=0;i<n;++i){double x=i/double(n-1),value=(target==std::string("a")?1.:2.)+(target==std::string("a")?1.:.5)*x*(1-x);error=std::max(error,std::abs(field->read(mesh.getIdx(i+1,1,1),0)-value));}}require(error<1e-9,"External resource block numerical error");require(original->values()[mesh.getIdx(9,1,1)]==.2+2*1.25-.4*2.125,"Borrowed input was overwritten");std::cout<<"borrowed model STATE error="<<error<<" owned arrays="<<resources.arrays.size()<<'\n';}
 require(state.time==3 && state.step==7 && state.dt==.125,"Steady block changed physical clock");
}
std::vector<double> producerBoundary(bool current){
 const int n=17;const double h=1./(n-1),dt=.0001,end=.01,pi=std::acos(-1.);
 BuildRequest r;r.includeDefaultFluidPreset=false;r.composition.stateDeclared=true;r.composition.solutionVariables={"a","b"};SystemContribution c;
 auto a=symbol("a",unit),b=symbol("b",diff),D=symbol("D",diff);D.role=StateRole::Derived;D.evaluation=StateEvaluation::Lazy;D.initializationRequired=D.boundaryRequired=D.outputEligible=false;D.dependencies={"b"};auto producer=std::make_shared<StateProducer>();producer->id="test.boundary.twice";producer->expression=X::constantValue(2)*b.expression();producer->versions={StateViewKind::Physical,StateViewKind::Stage};auto reads=std::make_shared<int>(0);
 producer->bind=[reads](Field& geometry,const auto& inputs,State::DistributedFieldView& output)->std::function<void()>{require(inputs.size()==1,"Producer dependency missing");output.read=[&geometry,input=inputs[0],reads](int cell,int){int i,j,k;geometry.getIJK(cell,i,j,k);require(i>=geometry.NG() && i<geometry.NG()+geometry.NX() && j==geometry.NG() && k==geometry.NG(),"Coefficient producer accessed unowned ghost");++*reads;return 2*input.read(cell,0);};return {};};D.producer=producer;c.addState(a);c.addState(b);c.addState(D);
 c.addEquation({"A",X::op(Math::Op::Ddt,{a.expression()}),X::op(Math::Op::Diffusion,{D.expression(),a.expression()})});c.addEquation({"B",X::op(Math::Op::Ddt,{b.expression()}),X::constantValue(-.7,rate)*b.expression()});
 for(const auto& pair:{std::pair{"A","a"},std::pair{"B","b"}}){ExecutionScope call;call.kind=ExecutionKind::EquationCall;call.id=pair.first;call.step={pair.first,{pair.second},call.id};c.addExecution(call);}
 NumericalBinding method{"A","ExplicitCentral2Arithmetic"};c.bindInputs({"A",{{"D",current?StateVersion::Current:StateVersion::Stage}}});c.bindNumerics(method);c.bindNumerics({"B","ScalarDiffusionCentral2"});ExecutionScope commit;commit.kind=ExecutionKind::Commit;commit.id="commit";commit.order=100;c.addExecution(commit);r.userContributions.push_back(std::move(c));
 FDM::SolverConfig config;config.numerics.maxDeltaT=dt;config.numerics.timeRecipe=FDM::builtInTimeRecipe(FDM::TimeRecipeId::ClassicalRK4);auto system=System::build(config,r);System::validate(system);
 Field mesh;mesh.setupGeometry(n,1,1,1);for(int k=0;k<mesh.MZ();++k)for(int j=0;j<mesh.MY();++j)for(int i=0;i<mesh.MX();++i)mesh.X(i,j,k)=(i-1)*h;mesh.getAllSets()["Left"]={mesh.getIdx(1,1,1)};mesh.getAllSets()["Right"]={mesh.getIdx(n,1,1)};
 CaseConfig input;auto description=std::make_shared<Model::Description>();input.modelDescription=description;for(const auto& name:{"a","b"}){Model::FieldDescriptor f;f.name=f.symbol=name;f.initial=0.;if(name==std::string("a"))for(const auto& side:{"Left","Right"})f.boundaries[side]={{"type","zeroGradient"}};description->fields.push_back(f);}
 State::StateBundle state;state.patches={&mesh};auto resources=Application::Execution::Detail::bindScalarResources(mesh,state,system,input);require(resources.arrays.size()==2,"Producer gained a physical copy");auto av=state.distributed.select("a",State::HaloSyncStage::None).at(0),bv=state.distributed.select("b",State::HaloSyncStage::None).at(0);std::vector<double> initial(n),coefficient(n);for(int i=0;i<n;++i){double x=i*h;initial[i]=1+.2*std::cos(pi*x);coefficient[i]=.03+.01*x+.002*std::cos(pi*x);av->write(mesh.getIdx(i+1,1,1),0,initial[i]);bv->write(mesh.getIdx(i+1,1,1),0,coefficient[i]);}
 FDM::SolverServices services;services.scalarInstances=resources.instances;Time::RunControl control;control.endTime=end;control.writeByStep=true;control.writeIntervalSteps=1000;Application::Execution::executeEquations(config,system,state,services,control);require(state.step==100 && std::abs(state.time-end)<1e-14 && *reads>0,"Producer/clock not executed");
 double effective=2*(1-std::exp(-.7*end))/.7;if(current){effective=0;double z=-.7*dt,R=1+z+z*z/2+z*z*z/6+z*z*z*z/24,scale=1;for(int step=0;step<state.step;++step){effective+=2*dt*scale;scale*=R;}}
 // Independent matrix-exponential reference of the existing mirror stencil; no production AST/kernel is used.
 const auto apply=[&](const std::vector<double>& v){std::vector<double> out(n);for(int i=0;i<n;++i)for(int sign:{-1,1}){int j=i+sign;if(j<0)j=1;else if(j>=n)j=n-2;out[i]+=.5*(coefficient[i]+coefficient[j])*(v[j]-v[i])/(h*h);}return out;};
 auto reference=initial,term=initial;for(int order=1;order<80;++order){term=apply(term);double largest=0;for(int i=0;i<n;++i){term[i]*=effective/order;reference[i]+=term[i];largest=std::max(largest,std::abs(term[i]));}if(largest<1e-19)break;require(order<79,"Independent exponential series failed");}
 std::vector<double> actual(n);double error=0,bError=0,budget=0,initialBudget=0;for(int i=0;i<n;++i){actual[i]=av->read(mesh.getIdx(i+1,1,1),0);error=std::max(error,std::abs(actual[i]-reference[i]));bError=std::max(bError,std::abs(bv->read(mesh.getIdx(i+1,1,1),0)-coefficient[i]*std::exp(-.7*end)));double weight=h*(i==0 || i==n-1?.5:1);budget+=weight*actual[i];initialBudget+=weight*initial[i];}require(error<2e-11 && bError<1e-14 && std::abs(budget-initialBudget)<1e-14,"Natural-boundary producer coefficient reference/budget failed");
 const auto name=current?"producerCurrentBoundary":"producerStageBoundary";auto folder=std::filesystem::path(__FILE__).parent_path()/"boundary/cases"/name/"result";std::filesystem::create_directories(folder);std::ofstream output(folder/"metrics.json");output<<std::setprecision(17)<<"{\"error\":"<<error<<",\"budgetError\":"<<std::abs(budget-initialBudget)<<",\"time\":"<<state.time<<",\"steps\":"<<state.step<<",\"coefficientPhysicalReads\":"<<*reads<<",\"arrays\":"<<resources.arrays.size()<<"}";std::cout<<name<<" error="<<error<<" budget="<<budget-initialBudget<<" physical reads="<<*reads<<'\n';return actual;
}

}
int main(int argc,char** argv){try{Parallel::Backend::MPIBackend mpi(argc,argv,true);run(false,false);run(true,false);run(false,true);auto stage=producerBoundary(false),current=producerBoundary(true);double difference=0;for(std::size_t i=0;i<stage.size();++i)difference=std::max(difference,std::abs(stage[i]-current[i]));require(difference>1e-9,"Current coefficient silently consumed Stage data");FDM::SolverConfig config;rejects([&]{System::build(config,request(true));},"layout");rejects([&]{System::build(config,request(false,true));},"version");std::cout<<"A8 real STATE borrowing, lifetime, whole solve and transactional failure passed\n";return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
