#define main standaloneReferenceMain
#include "test/test_standaloneScalar.cpp"
#undef main
#include "app/application/adapters/SF_scalarResources.h"
#include "infrastructure/io/SF_caseConfig.h"
namespace {
System::ResolvedSimulationSystem equations(FDM::SolverConfig& config) {
    using namespace System;using E=FormulaExpr;
    auto input=request("tracer",0.,true);
    auto& c=input.userContributions.front();
    const auto rate=Math::signature(0,{0,0,-1,0,0,0,0});
    c.registeredEquations.front()={"transport.tracer",E::op("ddt",{E::symbol("tracer",Math::signature(0,{}))}),E::op("source",{E::symbol("Sc",rate)})};
    StateSymbol second;second.id=second.name=second.storageKey="integral";second.signature=Math::signature(0,{});second.storageBinding=StorageBinding::NamedDistributed;second.restartEligible=false;c.addState(SF::Test::state(second));
    c.addEquation({"integrate",E::op("ddt",{E::symbol("integral",second.signature)}),E::multiply(E::constantValue(1,rate),E::symbol("tracer",second.signature))});
    c.execution.pop_back();ExecutionScope step;step.kind=ExecutionKind::EquationCall;step.id="second";step.step={"integrate",{"integral",TargetKind::Physical},"second"};c.addExecution(step);ExecutionScope commit;commit.kind=ExecutionKind::Commit;commit.id="commit";commit.order=1000000;c.addExecution(commit);c.bindNumerics({"integrate","ScalarDiffusionCentral2"});input.composition.solutionVariables.push_back("integral");
    auto result=build(config,input);
    // A port identity is a compiled resource reference, not an occurrence spelling.
    result.solvePlan.compiledProgram.steps.front().physicalResources.front().id="source.port";
    return result;
}
void sourceBinding() {
    FDM::SolverConfig config;config.numerics.maxDeltaT=.025;config.numerics.timeRecipe=FDM::builtInTimeRecipe(FDM::TimeRecipeId::ClassicalRK4);
    auto system=equations(config);Field mesh;mesh.setupGeometry(5,1,1,1);
    for(int k=0;k<mesh.MZ();++k)for(int j=0;j<mesh.MY();++j)for(int i=0;i<mesh.MX();++i)mesh.X(i,j,k)=(i-1)*.25;
    auto model=std::make_shared<Model::Description>();for(auto name:{"tracer","integral"}){Model::FieldDescriptor d;d.name=d.symbol=name;d.initial=0.;model->fields.push_back(d);}CaseConfig input;input.modelDescription=model;
    State::StateBundle state;state.patches={&mesh};
    std::vector<double> times;auto port=std::make_shared<FDM::ScalarEquationData>();port->target="tracer";port->sources["Sc"]=[&](double,double,double,double t){times.push_back(t);return t;};
    const auto requirement=system.solvePlan.compiledProgram.steps.front().physicalResources.front();
    auto borrowed=System::ResourceView::bind(requirement,port,*port,"native.scalar");
    auto resources=Application::Execution::Detail::bindScalarResources(mesh,state,system,input,0,{borrowed});
    require(resources.arrays.size()==2 && resources.ports.size()==1 && resources.instances.size()==2,"Compiled needs did not reuse source port");
    FDM::SolverServices services;services.scalarInstances=resources.instances;
    Time::RunControl control;control.endTime=.2;control.writeByStep=true;control.writeIntervalSteps=100;
    Application::Execution::executeEquations(config,system,state,services,control);
    const auto a=state.distributed.select("tracer",State::HaloSyncStage::None).front();const auto b=state.distributed.select("integral",State::HaloSyncStage::None).front();
    for(int i=1;i<=5;++i){int cell=mesh.getIdx(i,1,1);require(std::abs(a->read(cell,0)-.02)<1e-15,"Source did not use stage time");require(std::abs(b->read(cell,0)-.2*.2*.2/6)<1e-15,"Cross-equation synchronous input wrong");}
    require(state.step==8 && std::abs(state.time-.2)<1e-14,"Binder added a clock or timestep");
    require(times.size()==160 && times.at(0)==0 && std::abs(times.at(5)-.0125)<1e-15 && std::abs(times.at(15)-.025)<1e-15,"Source cadence/stage time wrong");
    const auto before=a->read(mesh.getIdx(2,1,1),0);
    auto reused=Application::Execution::Detail::bindScalarResources(mesh,state,system,input,0,{borrowed});require(reused.arrays.empty() && a->read(mesh.getIdx(2,1,1),0)==before,"Rebind copied/reset physical storage");
    // Failure must leave no dangling STATE registration, even after temporary allocation.
    for(int mode=0;mode<5;++mode){State::StateBundle bad;bad.patches={&mesh};auto p=std::make_shared<FDM::ScalarEquationData>(*port);auto r=requirement;if(mode==0)p->target="integral";if(mode==1)p->sources.clear();if(mode==2)r.writable=true;if(mode==3)r.owner="another.owner";
        auto view=System::ResourceView::bind(r,p,*p,r.owner);std::vector<System::ResourceView> views{view};if(mode==4)views.push_back(view);
        rejects([&]{Application::Execution::Detail::bindScalarResources(mesh,bad,system,input,0,views);});require(bad.distributed.select("tracer",State::HaloSyncStage::None).empty(),"Failed bind published dangling storage");}
    // A registered resource must serve the required identity, not just have its type.
    State::StateBundle missing;missing.patches={&mesh};auto wrong=requirement;wrong.id="unrelated";auto unrelated=System::ResourceView::bind(wrong,port,*port,wrong.owner);rejects([&]{Application::Execution::Detail::bindScalarResources(mesh,missing,system,input,0,{unrelated});});
    for(int mode=0;mode<3;++mode) {
        ScalarField backing;backing.setupLike(mesh,"tracer");auto view=State::scalarView("tracer",0,mesh,backing);
        if(mode==0)view.components=3;
        if(mode==1)view.write={};
        if(mode==2)view.requireReadable=[]{throw std::runtime_error("Backing lifetime expired");};
        State::StateBundle bad;bad.patches={&mesh};
        rejects([&]{bad.distributed.add(std::move(view));Application::Execution::Detail::bindScalarResources(mesh,bad,system,input,0,{borrowed});});
        require(bad.distributed.select("integral",State::HaloSyncStage::None).empty(),"Bad STATE shape/access/lifetime published another target");
    }
    // Runtime weak borrows cannot silently extend a user callback owner's lifetime.
    port.reset();times.clear();control.startTime=state.time;control.endTime=.225;
    bool expired=false;try{Application::Execution::executeEquations(config,system,state,services,control);}catch(const std::exception& error){expired=std::string(error.what()).find("Expired scalar resource owner: source.port")!=std::string::npos;}require(expired,"Expired source owner not accurately rejected");require(state.step==8 && std::abs(state.time-.2)<1e-14,"Expired resource committed a step");
    std::cout<<"source pair: 5x1x1 nodes, RK4 dt=.025 end=.2; a_t=t, b_t=a; a=.02 b=.0013333333333333333; eight physical commits\n";
}
}
int main(){try{sourceBinding();std::cout<<"Compiled field needs, borrowed source identity/lifetime, transactional binding and synchronous inputs passed\n";return 0;}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
