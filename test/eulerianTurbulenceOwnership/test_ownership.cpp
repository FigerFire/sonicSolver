#include "solver/algorithm/eulerian/equation/SF_equations.h"
#include "solver/boundary/SF_scalarBoundary.h"
#include "core/state/SF_distributedField.h"
#include "core/state/SF_stateBundle.h"
#include <algorithm>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace SF;
#include "referenceBoundary.h"
namespace {
void require(bool value,const char* message) {if(!value) throw std::runtime_error(message);}
template<class F> void rejects(F f) {bool failed=false;try{f();}catch(const std::exception&){failed=true;}require(failed,"Invalid ownership/epoch accepted.");}
struct Phases final:Physics::PhaseSystems::PhaseSystem {
    Field field;ScalarField pressure;std::vector<Physics::PhaseSystems::PhaseState> values;
    Physics::PhaseSystems::PhaseEquationSources source;
    Physics::Multiphase::MultiPhaseConfig configuration;
    Physics::Multiphase::PhaseProperties property;
    Phases() {
        field.setup(7,7,2,3,5);Math::resetActiveDirections();Math::deactivateDirection(2);
        for(int cell=0;cell<field.TotalSize();++cell) {
            int i,j,k;field.getIJK(cell,i,j,k);field.X(i,j,k)=i;field.Y(i,j,k)=j;field.Z(i,j,k)=k;
            field.Jac(i,j,k)=1;field.XiX(i,j,k)=1;field.EtY(i,j,k)=1;field.ZeZ(i,j,k)=1;
            if(i==3)field.getAllSets()["left"].push_back(cell);
            if(i==9)field.getAllSets()["right"].push_back(cell);
            if(j==3)field.getAllSets()["bottom"].push_back(cell);
            if(j==9)field.getAllSets()["top"].push_back(cell);
            if(k==3 || k==4)field.getAllSets()["empty"].push_back(cell);
        }
        pressure.setupLike(field,"p",1);property.viscosity=.001;
        for(const char* name:{"water","air"}) {
            values.emplace_back();auto& phase=values.back();phase.setupLike(field,name);
            std::fill(phase.primary.phaseMass.values().begin(),phase.primary.phaseMass.values().end(),values.size()*2);
            std::fill(phase.primitive.density.values().begin(),phase.primitive.density.values().end(),5);
            std::fill(phase.primitive.alpha.values().begin(),phase.primitive.alpha.values().end(),.4);
            std::fill(phase.primitive.temperature.values().begin(),phase.primitive.temperature.values().end(),300);
        }
    }
    void initialize(const Field&) override {} void recoverPrimitiveState() override {}
    void reconstructReferencePhaseMass() override {} void computeInterphase(double,const std::vector<Physics::PhaseSystems::PhaseVectorField>&) override {}
    void validateState(const std::string&) const override {}
    const Field& geometry() const override{return field;}
    ScalarField& sharedPressure() override{return pressure;}const ScalarField& sharedPressure() const override{return pressure;}
    std::vector<Physics::PhaseSystems::PhaseState>& phases() override{return values;}const std::vector<Physics::PhaseSystems::PhaseState>& phases() const override{return values;}
    Physics::PhaseSystems::PhaseEquationSources& sources() override{return source;}const Physics::PhaseSystems::PhaseEquationSources& sources() const override{return source;}
    size_t referencePhaseIndex() const override{return 0;}
    const Physics::Multiphase::PhaseProperties& phaseProperties(size_t) const override{return property;}
    const Physics::Multiphase::MultiPhaseConfig& config() const override{return configuration;}
};
FDM::TurbulenceConfig settings() {
    FDM::TurbulenceConfig c;c.enabled=true;c.family=FDM::TurbulenceFamily::RAS;c.model=FDM::TurbulenceModelKind::kEpsilon;c.phaseNames={"air","water"};
    c.scalars.kInitial={{"all",FIXED_VALUE,.01}};c.scalars.epsilonInitial={{"all",FIXED_VALUE,.02}};
    c.scalars.kBoundary={{"left",FIXED_VALUE,.01},{"right",ZERO_GRADIENT,0},{"bottom",ZERO_GRADIENT,0},{"top",ZERO_GRADIENT,0},{"empty",EMPTY,0}};
    c.scalars.epsilonBoundary=c.scalars.kBoundary;c.scalars.epsilonBoundary.front().value=.02;return c;
}
void boundaryTest(const Field& field) {
    ScalarField original;original.setupLike(field,"C");
    for(int cell=0;cell<field.TotalSize();++cell) original.values()[cell]=.3+cell*.0001;
    for(BCType law:{FIXED_VALUE,ZERO_GRADIENT,SYMMETRY,EMPTY}) {
        auto actual=original,reference=original;
        const std::vector<BCSetting<double>> data=law==EMPTY ? std::vector<BCSetting<double>>{{"empty",law,.7}}
            : std::vector<BCSetting<double>>{{"left",law,.7},{"right",law,.8},{"bottom",law,.9},{"top",law,1.0}};
        referenceBoundary(field,data,reference);Boundary::applyNodeScalarBoundary(field,data,actual);
        require(actual.values()==reference.values(),"Boundary/ghost differs from frozen node stencil.");
    }
    auto actual=original;
    rejects([&]{Boundary::applyNodeScalarBoundary(field,{{"left",static_cast<BCType>(999),0}},actual);});
    rejects([&]{Boundary::applyNodeScalarBoundary(field,{{"left",FIXED_VALUE,std::numeric_limits<double>::infinity()}},actual);});
}
}
int main() {
    Phases system;boundaryTest(system.field);auto config=settings();
    State::StateBundle multiple;Field other;multiple.patches={&system.field,&other};
    rejects([&]{(void)multiple.singlePatch();});
    Turbulence::EquationSystem model(config);model.initialize(system);
    EulerianEulerian::PhaseSolverWorkspace workspace;workspace.setupLike(system);
    FDM::PressureCorrectionConfig pressure;
    EulerianEulerian::PhaseEquationAssembler assembler(system,pressure,workspace);assembler.attachTurbulence(&model);
    require(workspace.turbulenceOldTime.size()==4,"Missing unique old-time workspace.");
    require(workspace.turbulenceOldTime[0].phaseIndex==1 && workspace.turbulenceOldTime[2].phaseIndex==0,"Selected order became phase storage order.");
    require(workspace.turbulenceOldTime[0].physical==&model.states()[0].kineticEnergy.variable,"Physical storage duplicated.");
    auto& target=model.states()[0].kineticEnergy.variable;
    auto view=State::scalarView("phase1.k",0,system.field,target);view.write(0,0,.2);require(target.values()[0]==.2,"STATE is not an alias.");
    const auto old=workspace.turbulenceOldTime[0].previousConserved.values();
    assembler.beginTurbulenceStep(5);
    rejects([&]{assembler.solveTurbulence(.01,5,0);});
    rejects([&]{assembler.publishTurbulenceCoefficients(6,0);});
    assembler.publishTurbulenceCoefficients(5,0);
    rejects([&]{assembler.solveTurbulence(.01,5,1);});
    auto& slot=workspace.turbulenceOldTime[0];const auto* physical=slot.physical;
    slot.physical=&slot.previousConserved;rejects([&]{assembler.solveTurbulence(.01,5,0);});slot.physical=physical;
    auto saved=std::move(slot.previousConserved);rejects([&]{assembler.solveTurbulence(.01,5,0);});slot.previousConserved=std::move(saved);
    slot.phaseIndex=0;rejects([&]{assembler.solveTurbulence(.01,5,0);});slot.phaseIndex=1;
    const int cell=system.field.getIdx(5,5,3);
    auto& sink=model.states()[0].kineticEnergy.implicitSink.values()[cell];const auto good=sink;
    sink=-1;rejects([&]{assembler.solveTurbulence(.01,5,0);});sink=good;
    auto& mass=system.values[1].primary.phaseMass.values()[cell];const auto goodMass=mass;
    mass=0;rejects([&]{assembler.solveTurbulence(.01,5,0);});mass=goodMass;
    auto& diffusivity=model.states()[0].kineticEnergy.diffusivity.values()[cell];const auto goodD=diffusivity;
    diffusivity=std::numeric_limits<double>::quiet_NaN();rejects([&]{assembler.solveTurbulence(.01,5,0);});diffusivity=goodD;
    auto& source=model.states()[0].kineticEnergy.explicitSource.values()[cell];const auto goodSource=source;
    source=std::numeric_limits<double>::quiet_NaN();
    bool linearRejected=false;
    try {assembler.solveTurbulence(.01,5,0);} catch(const std::exception& error) {
        linearRejected=std::string(error.what()).find("Eulerian turbulence solve failed")!=std::string::npos;
    }
    source=goodSource;require(linearRejected,"Invalid RHS did not reach the linear-session validation.");
    require(model.momentumEddyViscosity(1,cell)==model.states()[0].eddyViscosity.values()[cell]
        && model.energyEddyDiffusivity(1,cell)==model.states()[0].eddyViscosity.values()[cell]/config.turbulentPrandtl,
        "Momentum/energy feedback lost the original phase viscosity backing.");
    require(old==slot.previousConserved.values(),"Failed solve advanced committed old time.");
    rejects([&]{assembler.commitTimeLevel(5);});
    // Workspace commit uses the original physical mass, never selected-array ordinal or rho.
    workspace.commitTimeLevel(system);
    require(slot.previousConserved.values()[0]==system.values[1].primary.phaseMass.values()[0]*target.values()[0],"Wrong phase-mass old time.");
    Turbulence::EquationSystem bad([&]{auto c=config;c.coefficients.cMu=std::numeric_limits<double>::infinity();return c;}());bad.initialize(system);
    rejects([&]{bad.prepare(system);});
    EulerianEulerian::PhaseSolverWorkspace plainWorkspace;plainWorkspace.setupLike(system);
    EulerianEulerian::PhaseEquationAssembler plain(system,pressure,plainWorkspace);
    plain.beginTurbulenceStep(5);plain.commitTimeLevel(5);
    rejects([&]{plain.commitTimeLevel(5);});
    std::cout<<"Boundary all ghost layers, physical aliases, selected slots, old-time failure isolation and coefficient epochs Passed\n";
}
