#include "solver/algorithm/transport/SF_transportBinding.h"
#include "../sharedFluxTransport/fluid_fixture.h"
#include "solver/algorithm/SF_flowOperations.h"
#include "solver/algorithm/scalar/SF_scalarBinding.h"
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

    state.time=.3;state.step=17;
    auto realized=System::realizeState(system.executableSystem.state,system.runtime,state,system.solvePlan.compiledProgram.stateViews);
    const double limit=dt;
    auto flow=std::make_shared<SolverAlgorithm::FlowOperations>(config,system.executableSystem,system.numericalSystem,system.solvePlan,system.runtime,state,services,realized,limit);
    runtime.attachState(state);flow->prepare();Run::OpRegistry ops;flow->bindOperations(ops);
    auto f=flow->temporalParticipant();
    if(fault=="wrongProvidedSet") {
        auto altered=system.solvePlan;
        for(auto& call:altered.compiledProgram.steps)if(!call.coefficientPrepare.empty())call.coefficientStates.push_back("test.unimplemented-coefficient");
        for(auto& participant:altered.compiledProgram.temporalParticipants)if(!participant.coefficientPrepare.empty())participant.coefficientStates.push_back("test.unimplemented-coefficient");
        bool rejected=false;
        try{SolverAlgorithm::scalarParticipants(system.numericalSystem,altered,state,services,realized,limit,[flow](auto stamp){return flow->stageMassFlux(stamp);});}
        catch(const std::runtime_error& e){rejected=std::string(e.what()).find("workspace layout")!=std::string::npos;std::ofstream("negative.txt")<<e.what();}
        require(rejected,"A consistent compiled declaration overstated the actual coefficient workspace");
        return {};
    }
    auto scalars=SolverAlgorithm::scalarParticipants(system.numericalSystem,system.solvePlan,state,services,realized,limit,[flow](auto stamp){return flow->stageMassFlux(stamp);});
    std::ofstream trace("checkpoint-"+std::string(FDM::toString(id))+".txt");trace<<std::hexfloat;
    const auto physical=[&](const char* label) {
        trace<<label<<" step="<<state.step<<" time="<<state.time<<" dt="<<state.dt<<'\n';
        for(int i=g;i<g+n;++i) {
            trace<<i;for(int v=0;v<5;++v)trace<<' '<<mesh(i,g,g,v);
            trace<<' '<<turbulence.scalarFields().K(i,g,g)<<' '<<turbulence.scalarFields().Omega(i,g,g)<<' '<<turbulence.scalarFields().EddyMu(i,g,g)<<'\n';
        }
    };
    physical("initial");
    const auto flowRhs=f.rhs;
    f.rhs=[&](const auto& c) {
        flowRhs(c);const double t=state.time+system.numericalSystem.time.recipe.stage(c.stageIndex).abscissa*state.dt;
        trace<<"canonical-flux "<<state.step<<' '<<c.stageIndex<<' '<<t<<'\n';
        auto flux=flow->stageMassFlux({state.step,c.stageIndex,t});
        for(int i=g-1;i<g+n;++i)trace<<i<<' '<<flux.at(0,i,g,g)<<'\n';
        physical("flow-rhs");
    };
    std::vector<Time::TemporalCallbacks> callbacks{f};
    for(std::size_t slot=0;slot<scalars.size();++slot) {
        auto cb=scalars[slot];
        if(cb.coefficients) {
            auto original=cb.coefficients;cb.coefficients=[&,original](const auto& c) {
                original(c);trace<<"coefficients "<<state.step<<' '<<c.stageIndex<<' '<<state.time+system.numericalSystem.time.recipe.stage(c.stageIndex).abscissa*state.dt<<'\n';
                for(const auto* name:{"mu_t","SST.nu.k","SST.nu.omega","SST.production.k","SST.production.omega","SST.sink.k","SST.sink.omega","SST.cross"}) {
                    auto read=realized.stageReader(name,c.stageCount);trace<<name;
                    for(int i=g;i<g+n;++i)trace<<' '<<read(mesh.getIdx(i,g,g),0);trace<<'\n';
                }
            };
            auto publish=cb.publishCoefficients;cb.publishCoefficients=[&,publish]{publish();physical("final-coefficients");};
        }
        const auto rhs=cb.rhs;cb.rhs=[&,rhs,slot](const auto& c) {
            rhs(c);const auto name=slot==0?"k":"omega";
            trace<<"scalar-rhs "<<name<<' '<<state.step<<' '<<c.stageIndex<<'\n';
            const auto& workspace=realized.view(name,System::StateViewKind::Workspace).fields.front();
            auto primitive=realized.stageReader(name,c.stageCount);
            for(int i=g;i<g+n;++i){const int cell=mesh.getIdx(i,g,g);trace<<i<<' '<<workspace->read(cell,0)<<' '<<primitive(cell,0)*mesh(i,g,g,0)<<'\n';}
        };
        callbacks.push_back(std::move(cb));
    }
    Time::bindStageGroup(ops,system.solvePlan,system.numericalSystem,state,realized,limit,std::move(callbacks));
    for(const auto& b:system.runtime.operationBindings)if(b.operation==System::OpIds::TimeCommit)ops.bind(b.operation,b.provider,[&]{state.time+=state.dt;++state.step;});
    Run::PlanExecutor::validateBindings(system.solvePlan,ops);
    if(fault=="unpublished" || fault=="staleStage" || fault=="staleStep") {
        const auto& participants=system.solvePlan.compiledProgram.temporalParticipants;
        const auto invoke=[&](const std::string& op,const std::string& provider,int stage=0) {
            Run::ExecutionContext c;c.stageIndex=stage;c.stageCount=system.numericalSystem.time.recipe.stageCount();ops.invoke(op,c,provider);
        };
        const auto begin=[&] {
            for(const auto& p:participants)for(const auto& op:p.preparation)invoke(op,p.provider);
            invoke(System::TemporalOps::Dt,System::TemporalOps::Provider);
            for(const auto& p:participants)invoke(p.snapshot,p.provider);
        };
        const auto prepare=[&](int stage) {
            invoke(System::TemporalOps::Open,System::TemporalOps::Provider,stage);
            for(const auto& p:participants)invoke(p.prepareStage,p.provider,stage);
            invoke(System::TemporalOps::Ready,System::TemporalOps::Provider,stage);
        };
        const auto completeStage=[&](int stage) {
            prepare(stage);
            for(const auto& p:participants)if(!p.coefficientPrepare.empty())invoke(p.coefficientPrepare,p.provider,stage);
            for(const auto& p:participants)invoke(p.rhs,p.provider,stage);
            invoke(System::TemporalOps::RhsReady,System::TemporalOps::Provider,stage);
            for(const auto& p:participants)invoke(p.advance,p.provider,stage);
            invoke(System::TemporalOps::Close,System::TemporalOps::Provider,stage);
        };
        begin();
        if(fault=="staleStage")completeStage(0);
        if(fault=="staleStep") {
            for(int stage=0;stage<system.numericalSystem.time.recipe.stageCount();++stage)completeStage(stage);
            invoke(System::TemporalOps::PublishReady,System::TemporalOps::Provider);
            for(const auto& p:participants)invoke(p.publish,p.provider);
            for(const auto& p:participants)if(!p.coefficientPublish.empty())invoke(p.coefficientPublish,p.provider);
            for(const auto& b:system.runtime.operationBindings)if(b.operation==System::OpIds::TimeCommit)invoke(b.operation,b.provider);
            begin();
        }
        const int stage=fault=="staleStage"?1:0;prepare(stage);
        bool rejected=false;
        try{invoke(participants.front().rhs,participants.front().provider,stage);}
        catch(const std::runtime_error& e){rejected=std::string(e.what()).find("coefficient publication")!=std::string::npos;std::ofstream("negative.txt")<<e.what();}
        require(rejected,"Unpublished/previous stage/previous step coefficients were accepted");
        return {};
    }
    for(int step=0;step<3;++step) {
        const Run::PlanTraceContext c{state.step,state.time,&state.dt};Run::PlanExecutor::execute(system.solvePlan,ops,&c);physical("commit");
    }
    require(state.step==20 && std::abs(state.time-(.3+3*dt))<1e-14,"Coefficient harness changed physical clock");
    require(checked.publications==3*(system.numericalSystem.time.recipe.stageCount()+1),"Coefficient publication cadence changed");
    require(checked.viscousReads>0,"Mean-flow RHS did not read published stage viscosity");
    require(checked.samples==n*3*(system.numericalSystem.time.recipe.stageCount()+1),"Coefficient evaluation cadence changed");
    return {};
}
}
#ifndef SF_STAGE_COEFFICIENT_LIBRARY
int main(){for(auto id:{FDM::TimeRecipeId::ForwardEuler,FDM::TimeRecipeId::ClassicalRK4})sst(id,true,true);return 0;}

#endif
