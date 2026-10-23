#define main standaloneFixtureMain
#include "test/test_standaloneScalar.cpp"
#undef main
#include "app/application/adapters/SF_scalarResources.h"
#include "infrastructure/io/SF_caseConfig.h"
namespace {
using namespace System;using E=FormulaExpr;
BuildRequest modules(bool reverse=false,bool authored=true) {
    auto r=request("a",.03),b=request("b",.08);r.composition.solutionVariables.push_back("b");
    for(auto* c:{&r.userContributions.front(),&b.userContributions.front()}) {
        const auto id=c->states.front().id;c->records.clear();c->recordContribution("internal.diffusion","internal equation instance",id);
        c->execution.pop_back();
    }
    r.userContributions.push_back(b.userContributions.front());
    for(auto instance:{"first","second","reverse"}) {
        SystemContribution c;c.recordContribution("internal.reaction","additive local term",instance);
        const bool backwards=std::string(instance)=="reverse";
        c.extendMathematics({backwards?"transport.b":"transport.a"},E::multiply(
            E::constantValue(backwards?.3:.2,Math::signature(0,{0,0,-1,0,0,0,0})),E::symbol(backwards?"a":"b",Math::signature(0,{}))));
        r.modelContributions.push_back(c);
    }
    ExecutionScope commit;commit.kind=ExecutionKind::Commit;commit.id="commit";commit.order=100000;
    r.userContributions.back().execution.push_back(commit);
    if(authored) {
        ExecutionProgram how;how.root.kind=ExecutionKind::Sequence;how.root.id="declared.order";
        for(auto id:{"a","b"}) {ExecutionScope call;call.kind=ExecutionKind::EquationCall;call.step={"transport."+std::string(id),{id},"advance."+std::string(id)};how.root.children.push_back(call);}
        how.root.children.push_back(commit);r.authoredExecution=how;
        for(auto id:{"a","b"})r.authoredNumerics.push_back({"transport."+std::string(id),"ScalarDiffusionCentral2"});
    }
    if(reverse)std::reverse(r.userContributions.begin(),r.userContributions.end());
    return r;
}
std::vector<double> runModules(bool reverse) {
    FDM::SolverConfig config;config.numerics.maxDeltaT=.001;config.numerics.timeRecipe=FDM::builtInTimeRecipe(FDM::TimeRecipeId::ClassicalRK4);
    auto r=modules(reverse);auto system=build(config,r);validate(system);
    require(system.rawSystem.modifications.size()==3,"Identical intended additive terms lost");
    require(system.rawSystem.contributions.size()==5,"Module instance identity lost");
    for(const auto& c:system.rawSystem.contributions)require(!c.module.empty(),"Missing module category");
    Field geometry;geometry.setupGeometry(9,1,1,1);
    for(int k=0;k<geometry.MZ();++k)for(int j=0;j<geometry.MY();++j)for(int i=0;i<geometry.MX();++i)geometry.X(i,j,k)=(i-1)/8.;
    geometry.getAllSets()["Left"]={geometry.getIdx(1,1,1)};geometry.getAllSets()["Right"]={geometry.getIdx(9,1,1)};
    CaseConfig input;auto desc=std::make_shared<Model::Description>();input.modelDescription=desc;
    for(auto id:{"a","b"}) {Model::FieldDescriptor f;f.name=f.symbol=id;f.initial=std::string(id)=="a"?1.:2.;
        f.boundaries={{"Left",{{"type","zeroGradient"}}},{"Right",{{"type","zeroGradient"}}}};desc->fields.push_back(f);}
    State::StateBundle state;state.patches={&geometry};auto resources=Application::Execution::Detail::bindScalarResources(geometry,state,system,input);
    validatePhysicalResources(system.solvePlan.compiledProgram,resources.views);
    FDM::SolverServices services;services.scalarInstances=resources.instances;
    Time::RunControl control;control.endTime=.01;control.writeByStep=true;control.writeIntervalSteps=100;
    Application::Execution::executeEquations(config,system,state,services,control);
    require(state.step==10 && std::abs(state.time-.01)<1e-15,"New module changed clock");
    const double w=std::sqrt(.12),a=std::cosh(w*.01)+.8/w*std::sinh(w*.01),b=2*std::cosh(w*.01)+.3/w*std::sinh(w*.01);
    std::vector<double> result;
    for(const auto& id:{"a","b"}) {
        auto view=state.distributed.select(id,State::HaloSyncStage::None).at(0);
        for(int i=1;i<=9;++i) {double value=view->read(geometry.getIdx(i,1,1),0);require(std::abs(value-(id==std::string("a")?a:b))<1e-13,"Additive/cross-equation budget mismatch");result.push_back(value);}
    }
    std::ofstream("modules-"+std::to_string(reverse)+".txt")<<describe(system);
    return result;
}
std::string semanticSignature(bool reverse) {
    FDM::SolverConfig config;auto system=build(config,modules(reverse));
    std::map<std::string,std::string> math,state;
    for (const auto& e:system.rawSystem.registry.entries()) math[e.id]=canonicalFormula(e);
    for (const auto& v:system.rawSystem.state.symbols())state[v.id]=v.storageKey+":"+std::to_string(v.components)+":"+std::to_string(static_cast<int>(v.role));
    std::string result;
    for (const auto& [id,value]:math)result+=id+value;
    for (const auto& [id,value]:state)result+=id+value;
    for (const auto& binding:system.numericalSelection.bindings)result+=binding.equation+binding.method+binding.occurrence;
    const auto plan=[&](const auto& self,const SolvePlanNode& node)->void {
        result+=node.id+node.operation+node.provider+std::to_string(node.repetitions);
        for (const auto& child:node.children)self(self,child);
    };plan(plan,system.solvePlan.root);
    return result;
}
void conflicts() {
    FDM::SolverConfig config;
    auto r=modules(false,false);
    auto duplicate=r;duplicate.userContributions.push_back(duplicate.userContributions.front());rejects([&]{build(config,duplicate);});
    duplicate=r;auto c=duplicate.userContributions.front();c.records.clear();c.recordContribution("other","same symbol");duplicate.userContributions.push_back(c);rejects([&]{build(config,duplicate);});
    duplicate=r;duplicate.modelContributions[0].mathematicalExtensions[0].families={"missing"};rejects([&]{build(config,duplicate);});
    duplicate=r;duplicate.userContributions.front().bindNumerics({"transport.a","UnsupportedChoice"});rejects([&]{build(config,duplicate);});
    duplicate.authoredNumerics={{"transport.a","ScalarDiffusionCentral2"}};
    auto resolved=build(config,duplicate);require(std::any_of(resolved.numericalSelection.bindings.begin(),resolved.numericalSelection.bindings.end(),[](const auto& binding){return binding.equation=="transport.a" && binding.origin.kind==OriginKind::User;}),"Explicit choice lost provenance");
    duplicate.authoredNumerics.push_back(duplicate.authoredNumerics.front());rejects([&]{build(config,duplicate);});
    duplicate=r;duplicate.authoredNumerics={{"transport.a","UnknownMethod"}};rejects([&]{build(config,duplicate);});
    duplicate=r;
    for(int i=0;i<2;++i){SystemContribution c;c.recordContribution("closure","explicit closure owner",std::to_string(i));c.addClosure("c"+std::to_string(i),{},{},{"a"});duplicate.modelContributions.push_back(c);}
    rejects([&]{build(config,duplicate);});
    duplicate=modules();duplicate.authoredNumerics.clear();rejects([&]{build(config,duplicate);});
    // A producer is a single output declaration, not a competing implicit writer.
    duplicate=r;auto& state=duplicate.userContributions.front().states.front();state.shape=ValueShape::Vector;rejects([&]{build(config,duplicate);});
}
void boundResources() {
    CompiledExecutionProgram program;CompiledEquationCall call;call.physicalResources={{"port","owner",typeid(double),false}};program.steps.push_back(call);
    auto owner=std::make_shared<double>(3);auto view=ResourceView::bind(call.physicalResources.front(),owner,*owner,"owner");
    validatePhysicalResources(program,{view});rejects([&]{validatePhysicalResources(program,{});});rejects([&]{validatePhysicalResources(program,{view,view});});
    auto wrong=program;wrong.steps[0].physicalResources[0].owner="other";rejects([&]{validatePhysicalResources(wrong,{view});});
    owner.reset();rejects([&]{validatePhysicalResources(program,{view});});
}
}
int main() {
    try {require(semanticSignature(false)==semanticSignature(true),"Registration order changed canonical four-axis system/Plan");auto a=runModules(false),b=runModules(true);require(a==b,"Unrelated registration order changed authored HOW result");conflicts();boundResources();std::cout<<"Contribution composition and typed resource contracts passed\n";return 0;}
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
