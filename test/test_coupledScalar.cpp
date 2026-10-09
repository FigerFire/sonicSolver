#define main standaloneSingleMain
#include "test_standaloneScalar.cpp"
#undef main
#include "solver/algorithm/SF_singleFluidBinding.h"
#include "solver/algorithm/SF_singleFluidStepper.h"
namespace {
System::BuildRequest multiRequest(const std::string& a,const std::string& b,bool reverse=false) {
    auto r=request(a,.03,true),second=request(b,.08,true);
    r.userContributions.front().records.front().id+="."+a;
    second.userContributions.front().records.front().id+="."+b;
    r.userContributions.front().execution.pop_back();
    r.composition.solutionVariables.push_back(b);
    r.userContributions.push_back(second.userContributions.front());
    if (reverse) {
        std::reverse(r.userContributions.begin(),r.userContributions.end());
        r.userContributions.front().execution.front().order=-1;
    }
    return r;
}
System::BuildRequest coupledRequest(bool reverse=false,double gamma=.7) {
    auto r=multiRequest("a","b",reverse);using E=System::FormulaExpr;
    for(auto& c:r.userContributions) {
        auto& e=c.registeredEquations.front();const auto target=c.states.front().id;
        e.rhs=E::add(E::op("diffusion",{E::constantValue(.1),E::symbol(target)}),
            E::multiply(E::constantValue(gamma),E::symbol(target=="a"?"b":"a")));
    }
    return r;
}
std::array<Result,2> coupled(double dt,FDM::TimeRecipeId recipe,bool reverse=false,double start=0,double gamma=.7) {
    constexpr int n=9;const double h=1./8,pi=std::acos(-1.),end=start+.1;
    FDM::SolverConfig config;config.numerics.maxDeltaT=dt;config.numerics.timeRecipe=FDM::builtInTimeRecipe(recipe);
    auto system=System::build(config,coupledRequest(reverse,gamma));System::validate(system);
    Field mesh;mesh.setupGeometry(n,1,1,1);ScalarField a,b;a.setupLike(mesh,"a");b.setupLike(mesh,"b");
    for(int k=0;k<mesh.MZ();++k)for(int j=0;j<mesh.MY();++j)for(int i=0;i<mesh.MX();++i) {
        mesh.X(i,j,k)=(i-1)*h;a(i,j,k)=std::sin(pi*(i-1)*h);b(i,j,k)=0;
    }
    State::StateBundle state;state.patches={&mesh};state.time=start;state.step=start?17:0;
    state.distributed.add(State::scalarView("a",0,mesh,a));state.distributed.add(State::scalarView("b",0,mesh,b));
    auto da=std::make_shared<FDM::ScalarEquationData>(),db=std::make_shared<FDM::ScalarEquationData>();da->target="a";db->target="b";
    for(auto* data:{da.get(),db.get()}) {
        for(auto& bc:data->boundary)bc.kind=FDM::ScalarBoundaryKind::ZeroGradient;
        for(int side=0;side<2;++side) {data->boundary[side].kind=FDM::ScalarBoundaryKind::FixedValue;
            data->boundary[side].value=[](double,double,double,double){return 0;};}
    }
    FDM::SolverServices services;services.scalarInstances={{"advance.b",db},{"advance.a",da}};
    Time::RunControl control;control.startTime=start;control.endTime=end;control.writeByStep=true;control.writeIntervalSteps=1000;
    Application::Execution::executeEquations(config,system,state,services,control);
    require(std::abs(state.time-end)<1e-14 && state.step==int(std::round(.1/dt))+(start?17:0),"Coupled clock/restart wrong");
    std::array<Result,2> results;const double lambda=4*std::pow(std::sin(pi*h/2),2)/(h*h);
    for(int i=1;i<=n;++i) {
        const double spatial=std::sin(pi*(i-1)*h),t=.1;
        const std::array<double,2> discrete{std::exp(-.1*lambda*t)*std::cosh(gamma*t)*spatial,std::exp(-.1*lambda*t)*std::sinh(gamma*t)*spatial};
        const std::array<double,2> exact{std::exp(-.1*pi*pi*t)*std::cosh(gamma*t)*spatial,std::exp(-.1*pi*pi*t)*std::sinh(gamma*t)*spatial};
        for(int v=0;v<2;++v) {
            const double value=v?b(i,1,1):a(i,1,1);results[v].values.push_back(value);
            results[v].discreteError=std::max(results[v].discreteError,std::abs(value-discrete[v]));
            results[v].error=std::max(results[v].error,std::abs(value-exact[v]));
        }
    }
    std::cout<<"coupled dt="<<dt<<" semi-discrete="<<results[0].discreteError<<','<<results[1].discreteError
        <<" continuous="<<results[0].error<<','<<results[1].error<<'\n';return results;
}
void stageAndBindingGuards() {
    FDM::SolverConfig config;config.numerics.maxDeltaT=.01;config.numerics.timeRecipe=FDM::builtInTimeRecipe(FDM::TimeRecipeId::ClassicalRK4);
    auto system=System::build(config,coupledRequest());
    Field mesh;mesh.setupGeometry(9,1,1,1);ScalarField a,b;a.setupLike(mesh,"a",1);b.setupLike(mesh,"b",2);
    for(int k=0;k<mesh.MZ();++k)for(int j=0;j<mesh.MY();++j)for(int i=0;i<mesh.MX();++i)mesh.X(i,j,k)=(i-1)/8.;
    State::StateBundle state;state.patches={&mesh};state.distributed.add(State::scalarView("a",0,mesh,a));state.distributed.add(State::scalarView("b",0,mesh,b));
    auto da=std::make_shared<FDM::ScalarEquationData>(),db=std::make_shared<FDM::ScalarEquationData>();da->target="a";db->target="b";
    for(auto* d:{da.get(),db.get()})for(auto& bc:d->boundary)bc.kind=FDM::ScalarBoundaryKind::ZeroGradient;
    FDM::SolverServices services;services.scalarInstances={{"advance.a",da},{"advance.b",db}};
    auto realized=System::realizeState(system.executableSystem.state,system.runtime,state,system.solvePlan.compiledProgram.stateViews);
    Run::OpRegistry ops;double limit=.01;
    SolverAlgorithm::bindSingleFluidOperations(ops,config,system.executableSystem,system.numericalSystem,system.solvePlan,system.runtime,state,services,realized,limit);
    auto readA=realized.stageReader("a",4),readB=realized.stageReader("b",4);
    rejects([&]{readA(mesh.getIdx(4,1,1),0);});
    const auto& p=system.solvePlan.compiledProgram.temporalParticipants;
    ops.invoke(System::TemporalOps::Dt,{});for(const auto& item:p)ops.invoke(item.snapshot,{});
    Run::ExecutionContext context;context.stageCount=4;context.stageIndex=1;
    rejects([&]{ops.invoke(System::TemporalOps::Open,context);});context.stageIndex=0;
    ops.invoke(System::TemporalOps::Open,context);
    for(const auto& item:p)ops.invoke(item.prepareStage,context);
    ops.invoke(System::TemporalOps::Ready,context);
    const int cell=mesh.getIdx(4,1,1);
    require(readA(cell,0)==1 && readB(cell,0)==2,"Stage cross reads do not see snapshot");
    rejects([&]{realized.requireStageRead(1,0);});rejects([&]{realized.requireStageRead(0,.5);});
    ops.invoke(p[0].rhs,context);
    rejects([&]{ops.invoke(p[0].advance,context);});rejects([&]{ops.invoke(System::TemporalOps::RhsReady,context);});
    ops.invoke(p[1].rhs,context);require(readA(cell,0)==1 && readB(cell,0)==2,"RHS published stage too early");
    ops.invoke(System::TemporalOps::RhsReady,context);rejects([&]{readA(cell,0);});
    for(const auto& item:p)ops.invoke(item.advance,context);
    rejects([&]{ops.invoke(p[0].advance,context);});ops.invoke(System::TemporalOps::Close,context);
    require(a(4,1,1)==1 && b(4,1,1)==2 && state.time==0 && state.step==0,"Partial stage published physical state/clock");
    rejects([&]{ops.invoke(System::TemporalOps::PublishReady,{});});
    double expectedA=1+.005*.7*2,expectedB=2+.005*.7;
    for(int s=1;s<4;++s) {
        context.stageIndex=s;ops.invoke(System::TemporalOps::Open,context);
        for(const auto& item:p)ops.invoke(item.prepareStage,context);ops.invoke(System::TemporalOps::Ready,context);
        const double expectedTime=system.numericalSystem.time.recipe.stage(s).abscissa*.01;
        realized.requireStageRead(s,expectedTime);
        require(std::abs(readA(cell,0)-expectedA)<1e-14 && std::abs(readB(cell,0)-expectedB)<1e-14,"Cross-provider Stage values/index do not match RK recipe");
        const double aRhs=.7*expectedB,bRhs=.7*expectedA;
        expectedA=1+system.numericalSystem.time.recipe.stage(s).incrementWeight*.01*aRhs;
        expectedB=2+system.numericalSystem.time.recipe.stage(s).incrementWeight*.01*bRhs;
        for(const auto& item:p)ops.invoke(item.rhs,context);ops.invoke(System::TemporalOps::RhsReady,context);
        for(const auto& item:p)ops.invoke(item.advance,context);ops.invoke(System::TemporalOps::Close,context);
    }
    // Data lifetime error is detected for all participants before any physical publication.
    db.reset();rejects([&]{ops.invoke(System::TemporalOps::PublishReady,{});});
    require(a(4,1,1)==1 && b(4,1,1)==2 && state.time==0 && state.step==0,"Failed publication wrote half a solution");
    db=std::make_shared<FDM::ScalarEquationData>(*da);db->target="b";
    services.scalarInstances={{"advance.a",da},{"advance.b",db}};
    for(int mode=0;mode<6;++mode) {
        auto invalid=services;
        if(mode==0)invalid.scalarInstances[0].occurrence="wrong";
        if(mode==1)invalid.scalarInstances.push_back(invalid.scalarInstances.front());
        if(mode==2)invalid.scalarInstances.pop_back();
        auto wrong=std::make_shared<FDM::ScalarEquationData>(*da);
        if(mode==3){wrong->target="wrong";invalid.scalarInstances[0].data=wrong;}
        if(mode==4)wrong->sources.clear();
        auto invalidPlan=system.solvePlan;
        if(mode==4)invalidPlan.compiledProgram.steps[0].stateUses.back().version=System::StateVersion::Current;
        if(mode==5)invalidPlan.compiledProgram.steps[0].stateUses.pop_back();
        auto views=System::realizeState(system.executableSystem.state,system.runtime,state,system.solvePlan.compiledProgram.stateViews);
        Run::OpRegistry invalidOps;
        rejects([&]{SolverAlgorithm::bindSingleFluidOperations(invalidOps,config,system.executableSystem,system.numericalSystem,invalidPlan,system.runtime,state,invalid,views,limit);});
    }
    for(int mode=0;mode<4;++mode) {
        auto invalid=system.solvePlan;
        if(mode==0)invalid.root.children.push_back(invalid.root.children.back());
        const auto corrupt=[&](const auto& self,System::SolvePlanNode& node)->void {
            if(node.kind==System::PlanNodeKind::StageLoop) {
                if(mode==1)++node.repetitions;
                if(mode==2)node.children.erase(node.children.begin()+4);
                if(mode==3)std::swap(node.children[1],node.children[5]);
            }
            for(auto& child:node.children)self(self,child);
        };corrupt(corrupt,invalid.root);
        require(!System::validateTemporalPlan(invalid,4).empty(),"Malformed phase/Commit falsely valid");
        Run::OpRegistry invalidOps;auto views=System::realizeState(system.executableSystem.state,system.runtime,state,system.solvePlan.compiledProgram.stateViews);
        rejects([&]{SolverAlgorithm::bindSingleFluidOperations(invalidOps,config,system.executableSystem,system.numericalSystem,invalid,system.runtime,state,services,views,limit);});
        require(state.time==0 && state.step==0 && a(4,1,1)==1 && b(4,1,1)==2,"Invalid plan published partial physical/clock state");
    }
}
}
int main() {
    std::cout.precision(17);
    for(auto recipe:{FDM::TimeRecipeId::ForwardEuler,FDM::TimeRecipeId::ClassicalRK4}) {
        std::array<std::vector<double>,2> error;
        for(double dt:{.02,.01,.005}) {
            auto result=coupled(dt,recipe);for(int v=0;v<2;++v)error[v].push_back(result[v].discreteError);
        }
        const double threshold=recipe==FDM::TimeRecipeId::ForwardEuler?1.9:14;
        for(int v=0;v<2;++v) {
            std::cout<<"coupled ratios "<<error[v][0]/error[v][1]<<' '<<error[v][1]/error[v][2]<<'\n';
            require(error[v][0]/error[v][1]>threshold && error[v][1]/error[v][2]>threshold,"Coupled temporal order wrong");
        }
        auto uncoupled=coupled(.005,recipe,false,0,0);
        require(uncoupled[0].values==run("a",9,.1,.005,.1,recipe).values,"Zero coupling did not reduce to independent diffusion");
        for(double value:uncoupled[1].values)require(value==0,"Zero coupling created partner state");
        auto normal=coupled(.005,recipe),reversed=coupled(.005,recipe,true),restart=coupled(.005,recipe,false,.3);
        for(int v=0;v<2;++v) {
            require(normal[v].values==reversed[v].values,"Coupled HOW order changed solution");
            for(std::size_t i=0;i<normal[v].values.size();++i) require(std::abs(normal[v].values[i]-restart[v].values[i])<1e-14,"Restart clock truncation changed solution beyond roundoff");
        }
    }
    FDM::SolverConfig config;config.numerics.maxDeltaT=.01;
    auto absent=coupledRequest();absent.userContributions[1].states.clear();rejects([&]{System::build(config,absent);});
    auto unscheduled=coupledRequest();unscheduled.userContributions[1].execution.erase(unscheduled.userContributions[1].execution.begin());rejects([&]{System::build(config,unscheduled);});
    auto shape=coupledRequest();shape.userContributions[1].states[0].components=3;rejects([&]{System::build(config,shape);});
    auto nonlinear=coupledRequest();using E=System::FormulaExpr;
    nonlinear.userContributions[0].registeredEquations[0].rhs=E::add(nonlinear.userContributions[0].registeredEquations[0].rhs,E::multiply(E::symbol("a"),E::symbol("b")));
    rejects([&]{System::build(config,nonlinear);});
    auto duplicate=coupledRequest();duplicate.userContributions[1].execution.front()=duplicate.userContributions[0].execution.front();duplicate.userContributions[1].execution.front().step.occurrence="another";
    rejects([&]{System::build(config,duplicate);});
    auto commits=coupledRequest();commits.userContributions.back().execution.push_back(commits.userContributions.back().execution.back());
    rejects([&]{System::build(config,commits);});
    stageAndBindingGuards();
    std::cout<<"coupled scalar synchronous stages passed\n";return 0;
}
