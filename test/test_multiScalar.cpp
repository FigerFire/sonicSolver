// Reuse the standalone C++ client and analytical helpers, not another driver.
#define main standaloneSingleMain
#include "test_standaloneScalar.cpp"
#undef main
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
std::array<std::vector<double>,2> runMulti(const std::string& a,const std::string& b,
        FDM::TimeRecipeId recipe,bool reverse=false,bool varied=false) {
    constexpr int n=33;const double dt=.0005,end=.1,pi=std::acos(-1.);
    FDM::SolverConfig config;config.numerics.maxDeltaT=dt;
    config.numerics.timeRecipe=FDM::builtInTimeRecipe(recipe);
    auto system=System::build(config,multiRequest(a,b,reverse));System::validate(system);
    require(system.solvePlan.compiledProgram.temporalParticipants.size()==2,"Multi-instance owner count");
    Field mesh;mesh.setupGeometry(n,1,1,1);ScalarField first,second;
    first.setupLike(mesh,a);second.setupLike(mesh,b);
    for(int k=0;k<mesh.MZ();++k)for(int j=0;j<mesh.MY();++j)for(int i=0;i<mesh.MX();++i) {
        const double x=(i-1)/32.;mesh.X(i,j,k)=x;
        first(i,j,k)=std::sin(pi*x);second(i,j,k)=varied?2*std::cos(pi*x):std::sin(pi*x);
    }
    State::StateBundle state;state.patches={&mesh};
    state.distributed.add(State::scalarView(a,0,mesh,first));
    state.distributed.add(State::scalarView(b,0,mesh,second));
    auto da=std::make_shared<FDM::ScalarEquationData>(),db=std::make_shared<FDM::ScalarEquationData>();da->target=a;db->target=b;
    for(auto* data:{da.get(),db.get()}) {
        for(auto& bc:data->boundary) bc.kind=FDM::ScalarBoundaryKind::ZeroGradient;
        for(int side=0;side<2;++side) {data->boundary[side].kind=FDM::ScalarBoundaryKind::FixedValue;data->boundary[side].value=[](double,double,double,double){return 0;};}
    }
    da->sources["Sc"]=[pi](double x,double,double,double t){return (.03*pi*pi-1)*std::exp(-t)*std::sin(pi*x);};
    if(varied) {
        for(auto& bc:db->boundary)bc.kind=FDM::ScalarBoundaryKind::ZeroGradient;
        db->sources["Sc"]=[](double,double,double,double){return 0;};
    } else db->sources["Sc"]=[pi](double x,double,double,double t){return (.08*pi*pi-1)*std::exp(-t)*std::sin(pi*x);};
    FDM::SolverServices services;services.scalarInstances={{"advance."+b,db},{"advance."+a,da}};
    Time::RunControl control;control.endTime=end;control.writeByStep=true;control.writeIntervalSteps=100000;
    Application::Execution::executeEquations(config,system,state,services,control);
    require(mesh.NVar()==0 && !mesh.stateModel(),"Multi scalar inherited Q/EOS");
    require(state.step==200 && std::abs(state.time-end)<1e-14,"Repeated multi scalar clock commit");
    auto output=Application::Output::registeredStateFields(system.executableSystem.state,state,mesh);
    require(output.size()==2,"Multi scalar physical output missing");
    ResultWriterConfig settings;settings.caseDir="multi-scalar-output";settings.jobName="pair";
    std::filesystem::create_directories(settings.caseDir+"/result");
    ResultWriter(settings).save(mesh,state.time,output);
    std::ifstream written(settings.caseDir+"/result/pair_t0.1.vts");
    const std::string xml((std::istreambuf_iterator<char>(written)),{});
    require(xml.find("Name=\""+a+"\"")!=std::string::npos && xml.find("Name=\""+b+"\"")!=std::string::npos,"Multi scalar output file missing fields");
    std::ofstream(settings.caseDir+"/explain.txt")<<System::describe(system);
    std::array<std::vector<double>,2> result;
    for(int i=1;i<=n;++i) {result[0].push_back(first(i,1,1));result[1].push_back(second(i,1,1));}
    return result;
}
}
int main() {
    std::cout.precision(17);
    for(auto recipe:{FDM::TimeRecipeId::ForwardEuler,FDM::TimeRecipeId::ClassicalRK4}) {
        auto together=runMulti("a","b",recipe);
        auto singleA=run("a",33,.03,.0005,.1,recipe,true);
        auto singleB=run("b",33,.08,.0005,.1,recipe,true);
        require(together[0]==singleA.values && together[1]==singleB.values,"Independent multi scalar differs from separate execution");
        require(together==runMulti("a","b",recipe,true),"Declaration/HOW order changed independent solution");
        require(together==runMulti("arbitraryLeft","arbitraryRight",recipe),"Symbol renaming changed independent solution");
        const auto varied=runMulti("a","b",recipe,false,true);
        const auto differentB=run("b",33,.08,.0005,.1,recipe,false,{},true,2,0,false,true);
        require(varied[0]==singleA.values && varied[1]==differentB.values,"Independent initial/boundary/source ports mixed");
    }
    std::cout<<"multi-instance scalar Phase 1 passed\n";return 0;
}
