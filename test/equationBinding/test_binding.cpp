#include "solver/system/SF_lowering.h"
#include "solver/system/SF_systemBuilder.h"
#include "solver/system/SF_systemValidator.h"
#include "solver/system/SF_methodObjects.h"
#include "solver/system/SF_stateRealization.h"
#include <iostream>
using namespace SF;using namespace SF::System;
namespace {
void check(bool b,const char* why){if(!b)throw std::runtime_error(why);}
template<class F>void rejects(F f,const std::string& text){try{f();}catch(const std::exception& e){check(std::string(e.what()).find(text)!=std::string::npos,e.what());return;}throw std::runtime_error("Missing rejection "+text);}
BuildRequest make(){
 BuildRequest r;r.includeDefaultFluidPreset=false;r.composition.stateDeclared=true;r.composition.solutionVariables={"a","b","D"};SystemContribution c;
 for(auto id:{"a","b","D","D0"}){StateSymbol s;s.id=s.name=s.storageKey=id;s.signature=Math::signature(0,std::string(id)[0]=='D'?std::array<double,7>{0,2,-1,0,0,0,0}:std::array<double,7>{});s.storageBinding=StorageBinding::NamedDistributed;if(s.id=="D0"){s.constantValue=.001;s.runtimeStorageRequired=false;}c.addState(s);}
 auto a=c.states[0].expression(),b=c.states[1].expression(),D=c.states[2].expression(),D0=c.states[3].expression();using X=FormulaExpr;
 c.addEquation({"D",D,D0*(X::constantValue(1)+b)});c.addEquation({"G",X::op(Math::Op::Ddt,{a}),X::op(Math::Op::Diffusion,{D,a})});c.addEquation({"B",X::op(Math::Op::Ddt,{b}),X::constantValue(-.7,Math::signature(0,{0,0,-1,0,0,0,0}))*b});
 ExecutionScope stage;stage.kind=ExecutionKind::StageLoop;stage.id="stage";
 for(auto pair:{std::pair{"D","D"},std::pair{"G","a"},std::pair{"B","b"}}){ExecutionScope call;call.kind=ExecutionKind::EquationCall;call.id=std::string("call.")+pair.first;call.step={pair.first,{pair.second},call.id};stage.children.push_back(call);}
 c.addExecution(stage);ExecutionScope commit;commit.kind=ExecutionKind::Commit;commit.id="commit";commit.order=100000;c.addExecution(commit);
 c.bindNumerics({"D","LocalAffine"});c.bindNumerics({"G","ExplicitCentral2Arithmetic"});c.bindNumerics({"B","ScalarDiffusionCentral2"});r.userContributions.push_back(std::move(c));return r;
}
}
int main(){try{
 FDM::SolverConfig config;config.numerics.timeRecipe=FDM::builtInTimeRecipe(FDM::TimeRecipeId::ClassicalRK4);auto r=make();auto system=System::build(config,r);System::validate(system);check(Lowering::validateTemporalPlan(system.solvePlan,4).empty(),"Valid stage plan rejected");
 auto plan=system.solvePlan;auto stage=std::find_if(plan.root.children.begin(),plan.root.children.end(),[](const auto& n){return n.kind==PlanNodeKind::StageLoop;});auto action=std::find_if(stage->children.begin(),stage->children.end(),[](const auto& n){return n.occurrence=="call.D";});check(action!=stage->children.end(),"Missing declared stage action");stage->children.erase(action);check(!Lowering::validateTemporalPlan(plan,4).empty(),"Missing stage evaluation accepted");
 plan=system.solvePlan;stage=std::find_if(plan.root.children.begin(),plan.root.children.end(),[](const auto& n){return n.kind==PlanNodeKind::StageLoop;});action=std::find_if(stage->children.begin(),stage->children.end(),[](const auto& n){return n.occurrence=="call.D";});std::iter_swap(action,action+1);check(!Lowering::validateTemporalPlan(plan,4).empty(),"Wrong stage order accepted");
 plan=system.solvePlan;auto publication=std::find_if(plan.root.children.begin(),plan.root.children.end(),[](const auto& n){return n.occurrence=="call.D";});check(publication!=plan.root.children.end(),"No algebraic commit relation");plan.root.children.erase(publication);check(!Lowering::validateTemporalPlan(plan,4).empty(),"Missing physical relation publication accepted");
 auto readOnly=make();readOnly.userContributions[0].states[2].constantValue=.1;rejects([&]{System::build(config,readOnly);},"writable");
 auto wrong=make();std::swap(wrong.userContributions[0].execution[0].children[0],wrong.userContributions[0].execution[0].children[1]);rejects([&]{System::build(config,wrong);},"preceding Stage algebraic");
 auto block=make();ExecutionScope coupled;coupled.kind=ExecutionKind::BlockSolve;coupled.equations={"G","B","D"};coupled.unknowns={"a","b","D"};block.userContributions[0].execution={coupled};block.userContributions[0].numerics={{"*","ScalarDiffusionCentral2"}};rejects([&]{System::build(config,block);},"whole block requires an explicit Step binding");
 std::cout<<"Equation target, stage publication and block capability checks passed\n";return 0;
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
