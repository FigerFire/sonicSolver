#include "test/typedSymbols.h"
#include "solver/system/SF_systemBuilder.h"
#include "solver/system/SF_systemValidator.h"
#include "solver/system/SF_scalarMethod.h"
#include "solver/system/SF_systemPrinter.h"
#include "app/application/adapters/SF_scalarResources.h"
#include "app/application/execution/SF_flowLoop.h"
#include "app/application/output/SF_fields.h"
#include "infrastructure/io/SF_caseConfig.h"
#include "infrastructure/io/SF_resultWriter.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
using namespace SF;using namespace SF::System;using Expr=FormulaExpr;
namespace {
void require(bool v,const std::string& why){if(!v)throw std::runtime_error(why);}
template<class F> void rejects(F f,const std::string& expected){try{f();}catch(const std::exception& e){require(std::string(e.what()).find(expected)!=std::string::npos,e.what());std::cout<<"rejected "<<e.what()<<'\n';return;}throw std::runtime_error("Expected rejection: "+expected);}
const auto scalar=Math::signature(0,{}),vector=Math::signature(1,{0,1,-1,0,0,0,0});
Expr number(double v){return Expr::constantValue(v);}
Expr rate(double v){return Expr::constantValue(v,Math::signature(0,{0,0,-1,0,0,0,0}));}
Expr diffusion(double v,Expr x){return Expr::op(Math::Op::Diffusion,{Expr::constantValue(v,Math::signature(0,{0,2,-1,0,0,0,0})),x});}
BuildRequest request(const std::vector<std::string>& names,const std::vector<Expr>& rhs,bool vec=false,bool reverse=false){
 BuildRequest r;r.includeDefaultFluidPreset=false;r.composition.stateDeclared=true;r.composition.solutionVariables=names;SystemContribution c;c.recordContribution("test.explicit","A6 equations");
 for(std::size_t i=0;i<names.size();++i){StateSymbol s;s.id=s.name=s.storageKey=names[i];s.storageBinding=StorageBinding::NamedDistributed;s.signature=vec?vector:scalar;s.shape=vec?ValueShape::Vector:ValueShape::Scalar;s.components=vec?3:1;s.restartEligible=false;c.addState(s);
 c.addEquation({"equation."+s.id,Expr::op(Math::Op::Ddt,{s.expression()}),rhs[i]});ExecutionScope call;call.kind=ExecutionKind::EquationCall;call.id="advance."+s.id;call.step={"equation."+s.id,{s.id},call.id};c.addExecution(call);c.bindNumerics({"equation."+s.id,"ScalarDiffusionCentral2"});}
 if(reverse)std::reverse(c.execution.begin(),c.execution.end());ExecutionScope commit;commit.kind=ExecutionKind::Commit;commit.id="commit";commit.order=1000000;c.addExecution(commit);r.userContributions.push_back(std::move(c));return r;
}
struct Result{double error=0;std::vector<double> values;};
Result run(const std::string& id,BuildRequest req,int nodes,double dt,double end,FDM::TimeRecipeId time,
 std::function<double(int,int,double)> initial,std::function<double(int,int,double)> exact,bool stencil,bool vec=false,int ghost=1){
 FDM::SolverConfig config;config.numerics.maxDeltaT=dt;config.numerics.timeRecipe=FDM::builtInTimeRecipe(time);auto system=System::build(config,req);System::validate(system);
 require(system.runtime.report.status==RuntimeStatus::Runnable,"A6 plan not runnable");require(system.numericalSystem.requiredHaloWidth==(stencil?1:0),"Wrong compiled halo demand");
 Field mesh;mesh.setupGeometry(nodes,1,1,ghost);double h=1./(nodes-1);
 for(int k=0;k<mesh.MZ();++k)for(int j=0;j<mesh.MY();++j)for(int i=0;i<mesh.MX();++i)mesh.X(i,j,k)=(i-ghost)*h;
 mesh.getAllSets()["Left"]={mesh.getIdx(ghost,ghost,ghost)};mesh.getAllSets()["Right"]={mesh.getIdx(ghost+nodes-1,ghost,ghost)};
 CaseConfig input;auto description=std::make_shared<Model::Description>();input.modelDescription=description;
 for(const auto& name:req.composition.solutionVariables){Model::FieldDescriptor f;f.name=f.symbol=name;f.type=vec?Model::ValueType::Vector:Model::ValueType::Scalar;f.initial=vec?Model::Parameters::array({0.,0.,0.}):Model::Parameters(0.);
 if(stencil)for(const auto* side:{"Left","Right"})f.boundaries[side]={{"type","fixedValue"},{"value",f.initial}};description->fields.push_back(f);}
 State::StateBundle state;state.patches={&mesh};auto resources=Application::Execution::Detail::bindScalarResources(mesh,state,system,input);
 for(std::size_t q=0;q<req.composition.solutionVariables.size();++q){auto* view=state.distributed.select(req.composition.solutionVariables[q],State::HaloSyncStage::None).at(0);for(int i=0;i<nodes;++i)for(int c=0;c<(vec?3:1);++c)view->write(mesh.getIdx(i+ghost,ghost,ghost),c,initial(q,c,i*h));}
 const std::string path="cases/"+id;std::filesystem::create_directories(path+"/result");std::ofstream(path+"/explain.txt")<<System::describe(system);
 ResultWriterConfig output;output.caseDir=path;output.outputDir="result";output.jobName=id;ResultWriter writer(output);
 const auto fields=Application::Output::registeredStateFields(system.executableSystem.state,state,mesh);require(fields.front().components==(vec?3:1),"Vector output lost components");writer.save(mesh,0,fields);
 FDM::SolverServices services;services.scalarInstances=resources.instances;Time::RunControl control;control.endTime=end;control.writeByStep=true;control.writeIntervalSteps=100000;
 std::ofstream trace(path+"/clock.txt");trace.precision(17);Time::DriverCallbacks cb;cb.report=[&](int,const Time::AdvanceResult& result){trace<<result.time<<' '<<result.dt<<'\n';};
 Application::Execution::executeEquations(config,system,state,services,control,cb);writer.save(mesh,state.time,fields);
 require(std::abs(state.time-end)<1e-12 && mesh.NVar()==0,"Clock/Q changed");Result result;
 for(std::size_t q=0;q<req.composition.solutionVariables.size();++q){auto* view=state.distributed.select(req.composition.solutionVariables[q],State::HaloSyncStage::None).at(0);for(int i=0;i<nodes;++i)for(int c=0;c<(vec?3:1);++c){double v=view->read(mesh.getIdx(i+ghost,ghost,ghost),c);result.values.push_back(v);result.error=std::max(result.error,std::abs(v-exact(q,c,i*h)));}}
 std::ofstream metrics(path+"/metrics.json");metrics.precision(17);metrics<<"{\"error\":"<<result.error<<",\"time\":"<<state.time<<",\"steps\":"<<state.step<<",\"nodes\":"<<nodes<<",\"dt\":"<<dt<<",\"ghost\":"<<ghost<<",\"components\":"<<(vec?3:1)<<"}";
 std::cout<<id<<" error="<<result.error<<" time="<<state.time<<" steps="<<state.step<<'\n';return result;
}
}
int main(){try{
 std::cout.precision(17);const auto rk=FDM::TimeRecipeId::ClassicalRK4,fe=FDM::TimeRecipeId::ForwardEuler;const double pi=std::acos(-1.);
 auto a=Expr::symbol("a",scalar),b=Expr::symbol("b",scalar),c=Expr::symbol("c",scalar),v=Expr::symbol("arbitraryVector",vector);
 auto linear=request({"a","b"},{rate(.7)*b,rate(.3)*a});double z=std::sqrt(.21);auto init=[](int q,int,double){return q?2.:1.;};auto exact=[=](int q,int,double){return q?2*std::cosh(z)+.3/z*std::sinh(z):std::cosh(z)+1.4/z*std::sinh(z);};
 auto first=run("linearReaction",linear,5,.025,1,rk,init,exact,false,false,0);require(first.error<2e-9,"Linear reaction matrix exponential");
 auto reversed=request({"a","b"},{rate(.7)*b,rate(.3)*a},false,true);require(run("linearReverse",reversed,5,.025,1,rk,init,exact,false,false,0).values==first.values,"Cross-equation stage mixed old/new");
 auto logistic=request({"c"},{rate(2)*c*(number(1)-c)});
 for(auto method:{fe,rk}){std::vector<double> errors;for(double dt:{.1,.05,.025})errors.push_back(run("logistic"+std::string(method==fe?"Euler":"RK4")+std::to_string(dt),logistic,5,dt,1,method,[](int,int,double){return .2;},[](int,int,double){return 1./(1+4*std::exp(-2.));},false,false,0).error);require(errors[0]/errors[1]>(method==fe?1.7:12) && errors[1]/errors[2]>(method==fe?1.7:12),"Logistic temporal refinement");}
 for(int n:{17,33,65}){double h=1./(n-1),lambda=4*std::pow(std::sin(pi*h/2),2)/(h*h);double end=.05;auto multi=request({"a"},{diffusion(.03,a)+diffusion(.07,a)});auto s=run("multi"+std::to_string(n),multi,n,.0001,end,rk,[=](int,int,double x){return std::sin(pi*x);},[=](int,int,double x){return std::sin(pi*x)*std::exp(-.1*lambda*end);},true);require(s.error<2e-13,"Multiple diffusion semi-discrete reference");}
 {const int n=33;double h=1./(n-1),lambda=4*std::pow(std::sin(pi*h/2),2)/(h*h);auto duplicate=request({"a"},{diffusion(.05,a)+diffusion(.05,a)});
 require(run("duplicateDiffusion",duplicate,n,.0001,.05,rk,[=](int,int,double x){return std::sin(pi*x);},[=](int,int,double x){return std::sin(pi*x)*std::exp(-.1*lambda*.05);},true).error<2e-13,"Identical diffusion occurrences were deduplicated");}
 // Matrix [[.07,.02],[.02,.04]] has positive distinct eigenvalues. Independent analytic 2x2 exponential.
 auto cross=request({"a","b"},{diffusion(.07,a)+diffusion(.02,b),diffusion(.02,a)+diffusion(.04,b)});
 const int n=33;double h=1./(n-1),lambda=4*std::pow(std::sin(pi*h/2),2)/(h*h),end=.05;
 auto crossInitial=[=](int q,int,double x){return (q?2.:1.)*std::sin(pi*x);};
 auto crossExact=[=](int q,int,double x){double t=lambda*end,d=.025,mean=.055,sh=std::sinh(t*d)/d,ch=std::cosh(t*d);return std::sin(pi*x)*std::exp(-mean*t)*(q?2*ch-sh*(.02-.015*2):ch-sh*(.015+.02*2));};
 auto crossResult=run("crossDiffusion",cross,n,.0001,end,rk,crossInitial,crossExact,true);require(crossResult.error<3e-13,"Cross diffusion matrix exponential");
 auto signedCross=request({"a","b"},{diffusion(.07,a)+diffusion(-.02,b),diffusion(-.02,a)+diffusion(.04,b)});
 require(run("signedCross",signedCross,n,.0001,end,rk,crossInitial,[=](int q,int,double x){double sh=std::sinh(lambda*end*.025)/.025,ch=std::cosh(lambda*end*.025);return std::sin(pi*x)*std::exp(-.055*lambda*end)*(q?2*ch+sh*.05:ch+sh*.025);},true).error<3e-13,"Signed cross coefficients were altered");
 auto crossReverse=request({"a","b"},{diffusion(.07,a)+diffusion(.02,b),diffusion(.02,a)+diffusion(.04,b)},false,true);require(run("crossReverse",crossReverse,n,.0001,end,rk,crossInitial,crossExact,true).values==crossResult.values,"Cross diffusion stage order dependence");
 auto vectorEquation=request({"arbitraryVector"},{diffusion(.03,v)-rate(.3)*v},true);std::vector<double> space;
 for(int nodes:{17,33,65}){const double amplitude[3]={1,-2,.5};auto result=run("vector"+std::to_string(nodes),vectorEquation,nodes,.0001,.05,rk,[=](int,int c,double x){return amplitude[c]*std::sin((c+1)*pi*x);},[=](int,int c,double x){return amplitude[c]*std::sin((c+1)*pi*x)*std::exp(-(.03*(c+1)*(c+1)*pi*pi+.3)*.05);},true,true);space.push_back(result.error);}
 require(space[0]/space[1]>3.8 && space[1]/space[2]>3.8,"Vector spatial convergence");
 // A missing vector component or stencil input boundary fails before storage publication.
 {FDM::SolverConfig cfg;auto system=System::build(cfg,vectorEquation);Field mesh;mesh.setupGeometry(5,1,1,1);mesh.getAllSets()["Left"]={mesh.getIdx(1,1,1)};mesh.getAllSets()["Right"]={mesh.getIdx(5,1,1)};
 for(int mode=0;mode<3;++mode){CaseConfig input;auto desc=std::make_shared<Model::Description>();input.modelDescription=desc;Model::FieldDescriptor f;f.name=f.symbol="arbitraryVector";f.initial=Model::Parameters::array({1,2,3});for(auto side:{"Left","Right"})f.boundaries[side]={{"type","fixedValue"},{"value",f.initial}};
 if(mode==0)f.initial=Model::Parameters::array({1,2});if(mode==1)f.boundaries.erase("Right");if(mode!=2)desc->fields.push_back(f);State::StateBundle state;state.patches={&mesh};
 rejects([&]{Application::Execution::Detail::bindScalarResources(mesh,state,system,input);},mode==0?"explicit vector components":mode==1?"Missing scalar boundary":"Missing registered storage");require(state.distributed.select("arbitraryVector",State::HaloSyncStage::None).empty(),"Failed vector binder left a dangling owner");}}
 auto missing=request({"a"},{diffusion(.1,b)});rejects([&]{FDM::SolverConfig cfg;System::build(cfg,missing);},"Unknown base STATE");
 rejects([&]{auto wrong=diffusion(.1,a)+rate(.3)*v;(void)wrong;},"signature");
 rejects([&]{auto wrong=Expr::op(Math::Op::Diffusion,{rate(.1),a});Equation e{"bad",Expr::op(Math::Op::Ddt,{a}),wrong};},"signature");
 auto unsupported=request({"arbitraryVector"},{Expr::divide(Expr::op(Math::Op::Cross,{v,v}),Expr::constantValue(1,Math::signature(0,{0,1,0,0,0,0,0})))},true);
 rejects([&]{FDM::SolverConfig cfg;System::build(cfg,unsupported);},"no compatible numerical provider");
 std::cout<<"A6 explicit equation numerical closure passed\n";return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
