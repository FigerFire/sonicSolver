#include "solver/algorithm/transport/SF_transportBinding.h"
#include "fluid_fixture.h"
#include "models/turbulence/SF_turbulence.h"
#include "models/turbulence/SF_turbulenceSystemContribution.h"
#include "models/turbulence/RAS/SF_kOmegaSST.h"
#include "models/physics/equationRuntime/SF_mixtureEquationProvider.h"
#include "methods/numerics/viscous/SF_viscous.h"
#include "solver/discretization/diffusion/SF_formulaCentral2.h"
namespace {
struct ForbiddenImmersedBoundary : FDM::IImmersedBoundary {
    void apply(Field&,double,double) override{throw std::runtime_error("Unsupported stage IBM test boundary must never execute");}
};
struct CheckedTransport : FDM::ITransportModel {
    Turbulence::Manager& owner;mutable int samples=0,viscousReads=0;int publications=0;
    explicit CheckedTransport(Turbulence::Manager& m):owner(m){}
    void applyBoundary(const Field& f) override{owner.applyBoundary(f);}
    void refreshClosure(const Field&) override{throw std::runtime_error("Native SST entered legacy correct()");}
    double dynamicViscosity(const Field& f,int i,int j,int k,double mu)const override{++viscousReads;return owner.dynamicViscosity(f,i,j,k,mu);}
    std::vector<std::string> stageScalarTargets()const override{return owner.stageScalarTargets();}
    FDM::TransportStageCoefficients stageCoefficients(const Field& f,const FDM::TransportStageSamples& s,int i,int j,int k)const override{++samples;return owner.stageCoefficients(f,s,i,j,k);}
    FDM::ScalarEquationData stageScalarData(const Field& f,const std::string& n)const override{return owner.stageScalarData(f,n);}
    void publishStageViscosity(const Field& f,const std::vector<double>& v)override{++publications;owner.publishStageViscosity(f,v);}
};
std::vector<double> sst(FDM::TimeRecipeId id,bool active,bool diffusion=true,double dt=.0005,double end=.01,int manufactured=0,const std::string& fault={}) {
    const int n=33,g=3;const double h=1./(n-1),mu=.001;
    FDM::SolverConfig config;config.numerics.maxDeltaT=dt;config.numerics.cfl=2;
    config.numerics.timeRecipe=FDM::builtInTimeRecipe(id);config.numerics.viscousEnabled=diffusion;config.numerics.dynamicViscosity=mu;
    auto r=fluidRequest(false,false);r.singleFluidPreset=System::SingleFluidPresetSpec{diffusion};
    if(active){System::SystemContribution c;Turbulence::contribute(c,{"kOmegaSST",{},false,{},true});r.modelContributions.push_back(c);}
    if(manufactured)for(auto& equation:r.modelContributions.front().registeredEquations){using E=System::FormulaExpr;equation.rhs=E::op("diffusion",{E::multiply(E::symbol("rho"),E::constantValue(0)),E::symbol(equation.id)});}
    auto system=System::build(config,r);System::validate(system);
    require(std::find(system.runtime.report.requiredOperations.begin(),system.runtime.report.requiredOperations.end(),System::OpIds::TurbulenceAdvance)==system.runtime.report.requiredOperations.end(),"Native SST retained old advance authority");
    Field mesh;mesh.setup(n,1,1,g,5);auto model=Physics::FluidStateModel::makeSingleFluidPerfectGas(1.4,287.05,mu,.72);mesh.setStateModel(model);
    for(int cell=0;cell<mesh.TotalSize();++cell){int i,j,k;mesh.getIJK(cell,i,j,k);double x=(i-g)*h,u=manufactured?.3:.3+.05*std::sin(2*std::acos(-1.)*x),rho=manufactured==2?1+.1*std::sin(2*std::acos(-1.)*x):1.2;
        mesh.X(i,j,k)=x;mesh.Jac(i,j,k)=1/h;mesh.XiX(i,j,k)=1/h;mesh.EtY(i,j,k)=1;mesh.ZeZ(i,j,k)=1;mesh.wallDistance(i,j,k)=.1;
        mesh(i,j,k,0)=rho;mesh(i,j,k,1)=rho*u;mesh(i,j,k,4)=2.5+.5*rho*u*u;}
    mesh.getSetWritable("left")={mesh.getIdx(g,g,g)};mesh.getSetWritable("right")={mesh.getIdx(g+n-1,g,g)};
    for(int i=g;i<g+n;++i)mesh.getSetWritable("all").push_back(mesh.getIdx(i,g,g));
    FDM::TurbulenceConfig tc;tc.enabled=true;tc.family=FDM::TurbulenceFamily::RAS;tc.model=FDM::TurbulenceModelKind::kOmegaSST;tc.laminarDynamicViscosity=mu;
    tc.scalars.kInitial={{"all",FIXED_VALUE,.02}};tc.scalars.omegaInitial={{"all",FIXED_VALUE,2}};
    tc.scalars.kBoundary={{"left",ZERO_GRADIENT,0},{"right",ZERO_GRADIENT,0}};tc.scalars.omegaBoundary=tc.scalars.kBoundary;
    if(fault=="fixedBoundary"){for(auto& bc:tc.scalars.kBoundary){bc.type=FIXED_VALUE;bc.value=.02;}for(auto& bc:tc.scalars.omegaBoundary){bc.type=FIXED_VALUE;bc.value=2;}}
    Turbulence::Manager turbulence(tc);require(turbulence.initialize(mesh),"Model failed initialization");SolverAlgorithm::initializeTransport(mesh,turbulence);
    if(active)for(int i=g;i<g+n;++i){double x=(i-g)*h;turbulence.scalarFields().K(i,g,g)=manufactured==1?.02+.003*x:manufactured==2?.02:.02+.003*std::exp(-100*(x-.5)*(x-.5));turbulence.scalarFields().Omega(i,g,g)=manufactured==2?2:2+.1*x;}
    State::StateBundle state;state.patches={&mesh};state.stateModel=model;state.registerConservativeState();
    if(fault=="invalidStage")turbulence.scalarFields().K(g+5,g,g)=-1;
    if(active && fault!="missingStorage")turbulence.registerDistributed(state.distributed,0,mesh);
    CheckedTransport checked(turbulence);Equation::Coupling::CompositeTransportProvider composed;composed.setTurbulence(&checked);Execution::LocalRuntime runtime;HomogeneousExterior exterior;FDM::SolverServices services;services.executionRuntime=&runtime;services.boundaryPipeline=&exterior;if(active)services.transportModel=&composed;
    ForbiddenImmersedBoundary forbidden;if(fault=="immersed")services.immersed.boundary=&forbidden;
    Time::RunControl control;control.endTime=end;control.writeByStep=true;control.writeIntervalSteps=10000;
    if(active && !manufactured && fault.empty()) {
        FDM::TransportStageSamples samples={[](int,int,int){return .02;},[](int,int,int){return 2.;}};
        std::vector<double> muStage(mesh.TotalSize());for(int cell=0;cell<mesh.TotalSize();++cell){int i,j,k;mesh.getIJK(cell,i,j,k);muStage[cell]=.012;}
        checked.publishStageViscosity(mesh,muStage);std::vector<double> a(5),b(5);
        Viscous::TransportProperties molecular;molecular.mu=mu;auto turbulent=molecular;turbulent.transportModel=&checked;
        Viscous::viscousFluxAtFace(mesh,g+8,g,g,static_cast<Math::Dir>(0),Viscous::CentralOrder::SECOND,molecular,a);
        Viscous::viscousFluxAtFace(mesh,g+8,g,g,static_cast<Math::Dir>(0),Viscous::CentralOrder::SECOND,turbulent,b);
        double difference=0;for(int v=0;v<5;++v)difference=std::max(difference,std::abs(a[v]-b[v]));
        require(difference>1e-10,"Stage eddy viscosity did not alter actual viscous face flux");std::cout<<"viscous face flux feedback delta="<<difference<<'\n';checked.publications=0;
    }
    Application::Execution::executeEquations(config,system,state,services,control);
    require(std::abs(state.time-end)<1e-13 && state.step==int(std::round(end/dt)),"Native SST changed common clock");
    if(active){require(checked.publications==state.step*(system.numericalSystem.time.recipe.stageCount()+1),"Closure stage/final publication cadence wrong");require(!diffusion || checked.viscousReads>0,"NS did not consume stage viscosity");}
    std::vector<double> q;double minK=1e10,minW=1e10,maxMu=0;for(int i=g;i<g+n;++i){for(int v=0;v<5;++v)q.push_back(mesh(i,g,g,v));minK=std::min(minK,turbulence.scalarFields().K(i,g,g));minW=std::min(minW,turbulence.scalarFields().Omega(i,g,g));maxMu=std::max(maxMu,turbulence.scalarFields().EddyMu(i,g,g));}
    require(minK>0 && minW>0 && std::isfinite(maxMu),"SST physical fields invalid");
    if(fault=="fixedBoundary")for(int i:{g,g+n-1}){require(turbulence.scalarFields().K(i,g,g)==.02 && turbulence.scalarFields().Omega(i,g,g)==2,"Stage SST fixedValue boundary was not committed");}
    std::cout<<"sst "<<FDM::toString(id)<<" active="<<active<<" end="<<state.time<<" steps="<<state.step<<" publications="<<checked.publications<<" coefficientSamples="<<checked.samples<<" viscousReads="<<checked.viscousReads<<" kMin="<<minK<<" omegaMin="<<minW<<" muMax="<<maxMu<<'\n';
    if(active && !manufactured){std::ofstream out("sst-explain-"+std::string(FDM::toString(id))+".txt");out<<System::describe(system);}
    if(active)for(int i=g;i<g+n;++i){q.push_back(turbulence.scalarFields().K(i,g,g));q.push_back(turbulence.scalarFields().Omega(i,g,g));}
    return q;
}
void localKernel() {
    Field f;f.setup(5,1,1,3,5);f.setStateModel(Physics::FluidStateModel::makeSingleFluidPerfectGas(1.4,287.05,.001,.72));
    for(int cell=0;cell<f.TotalSize();++cell){int i,j,k;f.getIJK(cell,i,j,k);f(i,j,k,0)=1.2;f(i,j,k,4)=2.5;f.Jac(i,j,k)=1;f.XiX(i,j,k)=1;f.EtY(i,j,k)=1;f.ZeZ(i,j,k)=1;f.wallDistance(i,j,k)=.1;}
    FDM::TurbulenceConfig config;config.laminarDynamicViscosity=.001;
    FDM::TransportStageSamples samples={[](int,int,int){return .02;},[](int,int,int){return 2.;}};
    Turbulence::RAS::KOmegaSSTModel model;
    auto c=model.stageCoefficients(f,samples,config,4,3,3);
    Turbulence::ScalarFields old;old.resizeLike(f);for(int k=0;k<f.MZ();++k)for(int j=0;j<f.MY();++j)for(int i=0;i<f.MX();++i){old.K(i,j,k)=.02;old.Omega(i,j,k)=2;}
    double legacyMu=model.eddyDynamicViscosity(f,old,config,4,3,3);
    require(c.dynamicViscosity==legacyMu,"Stage viscosity changed legacy SST formula");
    auto settings=config;settings.enabled=true;settings.family=FDM::TurbulenceFamily::RAS;settings.model=FDM::TurbulenceModelKind::kOmegaSST;
    Turbulence::Manager once(settings);require(once.initialize(f),"Missing once-step fixture");SolverAlgorithm::initializeTransport(f,once);
    once.scalarFields()=old;
    const double dt=1e-5;SolverAlgorithm::updateOnceStepTransport(f,once,dt);
    const auto& updated=once.scalarFields();
    require(std::abs(updated.K(3,3,3)-(.02+dt*(c.production[0]-c.sink[0]*.02)/1.2))<1e-14,"SST k source math changed");
    require(std::abs(updated.Omega(3,3,3)-(2+dt*(c.production[1]-c.sink[1]*2+c.crossSource)/1.2))<1e-14,"SST omega source math changed");
    std::cout<<"SST zero-advection source/viscosity kernel preservation passed mu="<<legacyMu<<" F1="<<c.blending1<<" F2="<<c.blending2<<'\n';
}
void otherKernels() {
    Field f;f.setup(9,1,1,3,5);f.setStateModel(Physics::FluidStateModel::makeSingleFluidPerfectGas(1.4,287.05,.001,.72));
    for(int cell=0;cell<f.TotalSize();++cell){int i,j,k;f.getIJK(cell,i,j,k);double x=i-3;f(i,j,k,0)=1.2;f(i,j,k,1)=1.2*(.3+.05*x);f(i,j,k,4)=2.5+.5*1.2*std::pow(.3+.05*x,2);f.Jac(i,j,k)=1;f.XiX(i,j,k)=1;f.EtY(i,j,k)=1;f.ZeZ(i,j,k)=1;f.wallDistance(i,j,k)=10;}
    FDM::TurbulenceConfig config;config.laminarDynamicViscosity=.001;Turbulence::RAS::KOmegaSSTModel model;
    FDM::TransportStageSamples samples={[](int i,int,int){return .02+.001*(i-3);},[](int i,int,int){return 2+.1*(i-3);}};
    auto c=model.stageCoefficients(f,samples,config,7,3,3);double expected=2*(1-c.blending1)*1.2*.856*.001*.1/2.4;
    require(c.production[0]>0 && c.sink[0]>0 && c.crossSource>0,"SST nonzero production/sink/cross diffusion not exercised");require(std::abs(c.crossSource-expected)<1e-16,"SST cross diffusion sign/rho/F1 weight failed");
    System::FormulaValues v;v.spacing=[](int a){return a==0?1.:0.;};v.boundReads={[](int cell,int){return .02+.001*cell*cell;},[](int,int){return 1.2;},[](int,int){return 0.;},[=](int,int){return c.diffusivity[0];}};
    v.boundNeighbor=[&](std::size_t slot,int cell,int,int offset,int component){return v.boundReads.at(slot)(cell+offset,component);};
    double diffusion=Discretization::conservativeScalarCentral2(0,3,true).compileValue({})(v,4,0);
    require(std::abs(diffusion-.002*1.2*c.diffusivity[0])<1e-16,"SST concentration diffusion coefficient/sign failed");
    std::cout<<"SST production="<<c.production[0]<<" sink="<<c.sink[0]<<" cross="<<c.crossSource<<" F1="<<c.blending1<<" rhoWeightedDiffusion="<<diffusion<<'\n';
    auto bad=samples;bad[0]=[](int,int,int){return -1.;};rejects([&]{model.stageCoefficients(f,bad,config,7,3,3);});bad=samples;bad[1]=[](int,int,int){return std::numeric_limits<double>::quiet_NaN();};rejects([&]{model.stageCoefficients(f,bad,config,7,3,3);});
}
void negativeContracts() {
    FDM::SolverConfig config;config.numerics.maxDeltaT=.001;config.numerics.timeRecipe=FDM::builtInTimeRecipe(FDM::TimeRecipeId::ClassicalRK4);
    auto r=fluidRequest(false,false);System::SystemContribution c;Turbulence::contribute(c,{"kOmegaSST",{},false,{},true});r.modelContributions={c};
    const auto fail=[&](System::BuildRequest q){rejects([&]{auto s=System::build(config,q);System::validate(s);});};
    auto missing=r;missing.modelContributions[0].states.erase(missing.modelContributions[0].states.begin()+1);fail(missing);
    auto duplicate=r;duplicate.modelContributions[0].execution.push_back(c.execution[0]);fail(duplicate);
    auto invalid=r;invalid.modelContributions[0].registeredEquations[0].rhs=System::FormulaExpr::constantValue(0);fail(invalid);
    auto mismatch=r;mismatch.modelContributions[0].numerics[0].inputs={"k","epsilon"};fail(mismatch);
    auto mpi=r;mpi.parallel=true;fail(mpi);
    auto plan=System::build(config,r);const auto erase=[&](const auto& self,System::SolvePlanNode& node)->void{node.children.erase(std::remove_if(node.children.begin(),node.children.end(),[](const auto& c){return c.operation.find("scalar.stage.coefficients")!=std::string::npos;}),node.children.end());for(auto& child:node.children)self(self,child);};erase(erase,plan.solvePlan.root);require(!System::validateTemporalPlan(plan.solvePlan,4).empty(),"Plan accepted missing coefficient publication");
    rejects([&]{System::SystemContribution q;Turbulence::contribute(q,{"kEpsilon",{},false,{},true});});
    rejects([&]{sst(FDM::TimeRecipeId::ForwardEuler,true,true,.0005,.0005,0,"missingStorage");});
    rejects([&]{sst(FDM::TimeRecipeId::ForwardEuler,true,true,.0005,.0005,0,"invalidStage");});
    rejects([&]{sst(FDM::TimeRecipeId::ForwardEuler,true,true,.0005,.0005,0,"immersed");});
    (void)sst(FDM::TimeRecipeId::ClassicalRK4,true,true,.0005,.001,0,"fixedBoundary");
    auto pressure=r;pressure.singleFluidPreset.reset();pressure.composition.declared=true;pressure.composition.solutionVariables={"U","p"};pressure.composition.equations={"Continuity","Momentum"};pressure.composition.algorithm="PISO";pressure.composition.thermoDynamics.equationOfState="rhoConst";pressure.composition.thermoDynamics.transport="const";pressure.composition.thermoDynamics.constantDensity=1;
    fail(pressure);
    std::cout<<"native SST negative contracts passed\n";
}
void advectionAndDensity() {
    const double dt=.0005;
    for(int mode:{1,2}) {
        auto result=sst(FDM::TimeRecipeId::ForwardEuler,true,false,dt,dt,mode);double error=0;
        for(int cell=4;cell<29;++cell){double x=cell/32.;double k=mode==1?.02+.003*(x-.3*dt):.02,w=mode==1?2+.1*(x-.3*dt):2;
            error=std::max({error,std::abs(result[165+2*cell]-k),std::abs(result[165+2*cell+1]-w)});}
        require(error<1e-13,"Primitive SST targets did not use conservative shared advection/stage density");std::cout<<"SST targets manufactured mode="<<mode<<" advection/constant error="<<error<<'\n';
    }
}
void timeConvergence() {
    for(auto recipe:{FDM::TimeRecipeId::ForwardEuler,FDM::TimeRecipeId::ClassicalRK4}) {
        auto a=sst(recipe,true,true,.005,.04),b=sst(recipe,true,true,.0025,.04),reference=sst(recipe,true,true,.0003125,.04);
        double ea=0,eb=0;for(std::size_t i=0;i<a.size();++i){ea=std::max(ea,std::abs(a[i]-reference[i]));eb=std::max(eb,std::abs(b[i]-reference[i]));}
        require(ea>eb && eb>0,"Native SST time refinement did not reduce error");double order=std::log(ea/eb)/std::log(2.);std::cout<<"coupled SST time order="<<order<<" errors="<<ea<<","<<eb<<'\n';require(order>(recipe==FDM::TimeRecipeId::ForwardEuler?.8:3.5),"Coupled SST common stage time convergence failed");
    }
}

}
int main(){std::cout.precision(17);localKernel();otherKernels();negativeContracts();advectionAndDensity();timeConvergence();for(auto id:{FDM::TimeRecipeId::ForwardEuler,FDM::TimeRecipeId::ClassicalRK4}){auto laminar=sst(id,false),native=sst(id,true);double difference=0;for(std::size_t i=0;i<laminar.size();++i)difference=std::max(difference,std::abs(native[i]-laminar[i]));require(difference>1e-10,"SST feedback did not affect NS state");std::cout<<"NS feedback difference "<<difference<<'\n';}std::cout<<"stage SST execution passed\n";}
