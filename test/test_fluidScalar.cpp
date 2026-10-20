#define main standaloneSingleMain
#include "test_standaloneScalar.cpp"
#undef main
#include "models/physics/fluidStateModel/SF_factory.h"
#include "app/application/model/SF_output.h"
namespace {
struct HomogeneousExterior : FDM::IBoundaryPipeline {
    void prepare(const std::vector<Field*>& fields,double,double) override {
        for(auto* f:fields) {
            const int g=f->NG();
            for(int k=0;k<f->MZ();++k)for(int j=0;j<f->MY();++j)for(int i=0;i<f->MX();++i) {
                const int a=std::clamp(i,g,g+f->NX()-1),b=std::clamp(j,g,g+f->NY()-1),c=std::clamp(k,g,g+f->NZ()-1);
                if(i==a && j==b && k==c)continue;
                for(int v=0;v<f->NVar();++v)(*f)(i,j,k,v)=(*f)(a,b,c,v);
            }
            f->invalidateThermodynamicCache();
        }
    }
};
System::BuildRequest fluidRequest(bool scalar,bool feedback,bool reverse=false,double kappa=.5) {
    System::BuildRequest r;r.composition.stateDeclared=true;r.composition.solutionVariables={"rho","rhoU","rhoE"};
    r.singleFluidPreset=System::SingleFluidPresetSpec{false};
    if(!scalar)return r;
    auto extra=request("rhoC",.02,true).userContributions.front();using E=System::FormulaExpr;
    auto target=E::symbol("rhoC"),rho=E::symbol("rho"),c=E::divide(target,rho);
    extra.registeredEquations[0].lhs=E::add(E::op("ddt",{target}),E::op("div",{E::symbol("rhoU"),c}));
    extra.registeredEquations[0].rhs=E::add(E::op("diffusion",{E::multiply(rho,E::constantValue(.02)),c}),E::op("source",{E::symbol("Sc")}));
    extra.numerics[0].method="ScalarTransportUpwindCentral2";
    extra.execution.front().order=reverse?-1:100;
    extra.execution.pop_back();r.userContributions.push_back(extra);r.composition.solutionVariables.push_back("rhoC");
    if(feedback) {
        r.userContributions.back().extendMathematics({"momentum"},E::op("source",{E::symbol("linearScalarForce"),target,E::constantValue(kappa),E::constantValue(0)}));
        r.userContributions.back().extendMathematics({"energy"},E::op("source",{E::symbol("linearScalarWork"),target,E::constantValue(kappa),E::constantValue(0)}));
    }
    return r;
}
struct FluidResult {std::vector<double> q,scalar;double scalarError=0,momentumError=0,energyError=0;};
FluidResult fluid(double dt,FDM::TimeRecipeId recipe,bool scalar,bool feedback=false,bool reverse=false,bool nonlinearSource=false,double end=.1,double kappa=.5,bool output=false) {
    constexpr int n=17,g=3;const double h=1./(n-1),u0=.25,rho=1.2;
    FDM::SolverConfig config;config.numerics.maxDeltaT=dt;config.numerics.cfl=2;
    config.numerics.timeRecipe=FDM::builtInTimeRecipe(recipe);
    auto system=System::build(config,fluidRequest(scalar,feedback,reverse,kappa));System::validate(system);
    if(scalar) {
        const auto& participants=system.solvePlan.compiledProgram.temporalParticipants;
        require(participants.size()==2,"Flow fused owner split or scalar not in common group");
        require(participants.front().identity==(reverse?"advance.rhoC":"flow.conservative"),"Test did not reverse actual participant/RHS order");
    }
    Field mesh;mesh.setup(n,1,1,g,5);
    auto model=Physics::FluidStateModel::makeSingleFluidPerfectGas(1.4,287.05,0,.72);mesh.setStateModel(model);
    ScalarField tracer;tracer.setupLike(mesh,"rhoC");
    const double energy0=2.5+.5*rho*u0*u0;
    for(int cell=0;cell<mesh.TotalSize();++cell) {
        int i,j,k;mesh.getIJK(cell,i,j,k);mesh.X(i,j,k)=(i-g)*h;
        mesh.Jac(i,j,k)=1;mesh.XiX(i,j,k)=1/h;mesh.EtY(i,j,k)=1;mesh.ZeZ(i,j,k)=1;
        mesh(i,j,k,0)=rho;mesh(i,j,k,1)=rho*u0;mesh(i,j,k,4)=energy0;
        tracer(i,j,k)=feedback && kappa!=0?rho:rho*(1+mesh.X(i,j,k));
    }
    // Constant flow has exact zero divergence. Homogeneous exterior Q gives
    // a manufactured transport test without adding unrelated wall physics.
    State::StateBundle state;state.patches={&mesh};state.stateModel=model;state.registerConservativeState();
    auto data=std::make_shared<FDM::ScalarEquationData>();data->target="rhoC";
    for(auto& bc:data->boundary)bc.kind=FDM::ScalarBoundaryKind::ZeroGradient;
    if(!feedback || kappa==0)for(int side=0;side<2;++side) {
        data->boundary[side].kind=FDM::ScalarBoundaryKind::FixedValue;
        data->boundary[side].value=[=](double x,double,double,double t){return rho*(1+x-u0*t+(nonlinearSource?std::exp(t)-1-t:0));};
    }
    data->sources["Sc"]=[=](double,double,double,double t){return nonlinearSource?rho*(std::exp(t)-1):0;};
    Execution::LocalRuntime runtime;HomogeneousExterior exterior;FDM::SolverServices services;services.executionRuntime=&runtime;services.boundaryPipeline=&exterior;
    if(scalar){state.distributed.add(State::scalarView("rhoC",0,mesh,tracer));services.scalarInstances={{"advance.rhoC",data}};}
    Time::RunControl control;control.endTime=end;control.writeByStep=true;control.writeIntervalSteps=10000;
    Application::Execution::executeEquations(config,system,state,services,control);
    require(std::abs(state.time-end)<1e-13 && state.step==int(std::round(end/dt)),"Mixed providers changed common dt/clock");
    if(output) {
        System::StateRegistry outputSymbols;outputSymbols.add(system.executableSystem.state.at("rhoC"));
        auto fields=Application::Output::registeredStateFields(outputSymbols,state,mesh);
        require(std::any_of(fields.begin(),fields.end(),[](const auto& f){return f.name=="rhoC";}),"Flow + scalar physical output absent");
        ResultWriterConfig settings;settings.caseDir="fluid-scalar-output";settings.jobName="flowTracer";
        std::filesystem::create_directories(settings.caseDir+"/result");
        ResultWriter writer(settings);
        CaseConfig outputConfig;Application::ModelLoader::configureOutput(writer,outputConfig);
        writer.save(mesh,state.time,fields);
        std::ifstream written(settings.caseDir+"/result/flowTracer_t0.1.vts");
        const std::string xml((std::istreambuf_iterator<char>(written)),{});
        for(const auto* field:{"Density","Velocity","Pressure","rhoC"})
            require(xml.find(std::string("Name=\"")+field+"\"")!=std::string::npos,"Flow + scalar VTS missing physical array");
        std::ofstream("fluid-scalar-output/explain.txt")<<System::describe(system);
        require(System::describe(system).find("SYNCHRONOUS TEMPORAL GROUP")!=std::string::npos,"Explain omitted common stage group");
    }
    FluidResult result;
    for(int i=g;i<g+n;++i) {
        for(int v=0;v<5;++v)result.q.push_back(mesh(i,g,g,v));
        const double exactScalar=feedback && kappa!=0?rho*(nonlinearSource?std::exp(end)-end:1):rho*(1+(i-g)*h-u0*end+(nonlinearSource?std::exp(end)-1-end:0));
        result.scalar.push_back(tracer(i,g,g));result.scalarError=std::max(result.scalarError,std::abs(tracer(i,g,g)-exactScalar));
        const double velocity=u0+(feedback?kappa*(nonlinearSource?std::exp(end)-1-.5*end*end:end):0);
        result.momentumError=std::max(result.momentumError,std::abs(mesh(i,g,g,1)-rho*velocity));
        result.energyError=std::max(result.energyError,std::abs(mesh(i,g,g,4)-(energy0+.5*rho*(velocity*velocity-u0*u0))));
    }
    std::cout<<"fluid scalar="<<scalar<<" feedback="<<feedback<<" dt="<<dt<<" scalarError="<<result.scalarError
        <<" momentumError="<<result.momentumError<<" energyError="<<result.energyError<<'\n';return result;
}
}
int main() {
    std::cout.precision(17);
    for(auto recipe:{FDM::TimeRecipeId::ForwardEuler,FDM::TimeRecipeId::ClassicalRK4}) {
        auto reference=fluid(.005,recipe,false),passive=fluid(.005,recipe,true,false,false,false,.1,.5,true),reordered=fluid(.005,recipe,true,false,true,false,.1,.5,true);
        require(reference.q==passive.q && passive.q==reordered.q,"Passive scalar changed flow or group order");
        require(passive.scalar==reordered.scalar && passive.scalarError<1e-12,"Passive conservative advection did not use correct Stage flow");
        auto zeroFeedback=fluid(.005,recipe,true,true,false,false,.1,0);
        require(zeroFeedback.q==passive.q && zeroFeedback.scalar==passive.scalar,"Zero feedback did not reduce exactly to passive flow + scalar");
        auto feedback=fluid(.005,recipe,true,true),feedbackReverse=fluid(.005,recipe,true,true,true);
        require(feedback.q==feedbackReverse.q && feedback.scalar==feedbackReverse.scalar,"Feedback RHS ordering changed synchronous solution");
        require(feedback.momentumError<1e-12 && feedback.scalarError<1e-12,"Feedback momentum or conservative scalar incorrect");
        if(recipe==FDM::TimeRecipeId::ClassicalRK4)require(feedback.energyError<1e-12,"Mechanical work missed Stage velocity");
        std::vector<double> errors;
        for(double dt:{.02,.01,.005})errors.push_back(fluid(dt,recipe,true,false,false,true).scalarError);
        const double threshold=recipe==FDM::TimeRecipeId::ForwardEuler?1.9:14;
        require(errors[0]/errors[1]>threshold && errors[1]/errors[2]>threshold,"Multi-provider temporal convergence wrong");
        std::cout<<"multi-provider ratios "<<errors[0]/errors[1]<<' '<<errors[1]/errors[2]<<'\n';
        std::vector<double> momentum,energy;
        for(double dt:recipe==FDM::TimeRecipeId::ForwardEuler?std::vector<double>{.01,.005,.0025}:std::vector<double>{.04,.02,.01}) {
            auto value=fluid(dt,recipe,true,true,false,true,.4);
            momentum.push_back(value.momentumError);energy.push_back(value.energyError);
        }
        for(const auto* error:{&momentum,&energy}) {
            std::cout<<"feedback temporal ratios "<<(*error)[0]/(*error)[1]<<' '<<(*error)[1]/(*error)[2]<<'\n';
            require((*error)[0]/(*error)[1]>threshold && (*error)[1]/(*error)[2]>threshold,"Stage feedback momentum/work order wrong");
        }
    }
    std::cout<<"fluid + conservative scalar synchronous execution passed\n";
}
