#include "models/physics/phaseSystem/SF_eulerianEulerian.h"
#include "models/physics/phaseSystem/SF_phaseSource.h"
#include "models/physics/heat/SF_phaseWallHeat.h"
#include "models/physics/phaseChange/Lee/SF_Lee.h"
#include "models/physics/fluidStateModel/SF_factory.h"
#include "solver/system/SF_systemBuilder.h"
#include "solver/system/SF_pressureCoupling.h"
#include "solver/run/SF_planExecutor.h"
#include "core/system/SF_operationIds.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>

using namespace SF;
namespace PS = Physics::PhaseSystems;
namespace MP = Physics::Multiphase;
namespace {
void require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}
void near(double actual, double expected, const char* message) {
    if (!std::isfinite(actual) || !std::isfinite(expected)
        || std::abs(actual-expected)>2e-12*std::max({1.0,std::abs(actual),std::abs(expected)}))
        throw std::runtime_error(message);
}
MP::MultiPhaseConfig configuration(const std::string& model, bool enabled, double temperature) {
    MP::MultiPhaseConfig c;
    c.enabled=true; c.type="eulerianEulerian"; c.initialPressure=101325;
    c.eulerianEulerian.phaseNames={"liquid","vapor"}; c.eulerianEulerian.referencePhase="liquid";
    MP::PhasePairModelOptions pair;
    pair.name="liquid_vapor";pair.continuousPhase="liquid";pair.dispersedPhase="vapor";
    pair.particleDiameter=.001;pair.liftCoefficient=0;pair.virtualMassCoefficient=0;
    c.eulerianEulerian.phasePairs={pair};
    MP::PhaseProperties liquid, vapor;
    liquid.name="liquid"; liquid.role=MP::PhaseRole::Liquid; liquid.thermoModel="perfectLiquid";
    liquid.density=1000; liquid.volumeFraction=.8; liquid.viscosity=.001;
    liquid.specificHeat=liquid.Cv=4000;liquid.thermalConductivity=.6;
    vapor.name="vapor";vapor.role=MP::PhaseRole::Gas;vapor.thermoModel="perfectGas";
    vapor.density=c.initialPressure/(287.05*temperature);vapor.volumeFraction=.2;vapor.viscosity=1e-5;
    vapor.specificHeat=1000;vapor.Cv=713;vapor.gasConstant=287.05;vapor.thermalConductivity=.02;
    for (auto* p:{&liquid,&vapor}) {
        p->initialTemperature=temperature;p->initialTemperatureDeclared=true;
        p->initialVelocity=Vector3(.25,.1,0);
    }
    c.phases={liquid,vapor};c.alpha.phaseName="liquid";c.alpha.defaultValue=.8;
    auto& pc=c.phaseChange;
    pc.model=model;pc.enabled=enabled;pc.phaseNames={"liquid","vapor"};
    pc.latentHeat=2e6;pc.saturationTemperature=373;pc.evaporationCoefficient=.01;pc.condensationCoefficient=.02;
    pc.selections={{"departurediametermodel","explicit"},{"departurefrequencymodel","explicit"},
        {"nucleationsitedensitymodel","explicit"},{"convectiveheatfluxmodel","fixed"},
        {"quenchingheatfluxmodel","fixed"},{"walltemperaturemodel","explicit"}};
    const double vaporDensity=c.initialPressure/(vapor.gasConstant*temperature);
    const double bubbleVolume=std::acos(-1.0)*std::pow(.002,3)/6;
    pc.coefficients={{"departurediameter",.002},{"departurefrequency",20},
        {"nucleationsitedensity",.002/(bubbleVolume*vaporDensity*20)},
        {"convectiveheatflux",1000},{"quenchingheatflux",2000},
        {"maximumwallvaporfraction",.9},{"heatfluxbalancetolerance",1e-10}};
    return c;
}
Field geometry() {
    Field f;f.setup(7,7,2,3,5);
    f.setStateModel(Physics::FluidStateModel::makeSingleFluidPerfectGas(1.4,287.05,0,.72));
    Math::resetActiveDirections();Math::deactivateDirection(2);
    for (int cell=0;cell<f.TotalSize();++cell) {
        int i,j,k;f.getIJK(cell,i,j,k);
        f.X(i,j,k)=(i-3)*.25;f.Y(i,j,k)=(j-3)*.5;f.Z(i,j,k)=(k-3)*.1;
        f.Jac(i,j,k)=80;f.XiX(i,j,k)=4;f.EtY(i,j,k)=2;f.ZeZ(i,j,k)=10;
        f.CellFlag(i,j,k)=FLUID_CELL;
        f(i,j,k,0)=1;f(i,j,k,1)=f(i,j,k,2)=f(i,j,k,3)=0;f(i,j,k,4)=101325/.4;
        if(i>=3&&i<=9&&j>=3&&j<=9&&k>=3&&k<=4) {
            if(i==3)f.getAllSets()["left"].push_back(cell);
            if(i==9)f.getAllSets()["right"].push_back(cell);
            if(j==3)f.getAllSets()["bottom"].push_back(cell);
            if(j==9)f.getAllSets()["top"].push_back(cell);
        }
    }
    for(const auto* patch:{"left","right","bottom","top"}) f.markSolverBoundarySet(patch);
    return f;
}
using Snapshot = std::vector<std::vector<double>>;
Snapshot sourceSnapshot(const PS::PhaseEquationSources& s) {
    Snapshot result;
    for(size_t p=0;p<s.mass.size();++p) {
        result.push_back(s.mass[p].values());result.push_back(s.energy[p].values());
        for(int d=0;d<3;++d) result.push_back(s.momentum[p][d].values());
    }
    for(const auto* f:{&s.wallBoilingMass,&s.wallBoilingConvectiveHeat,&s.wallBoilingQuenchingHeat,&s.wallBoilingEvaporativeHeat})
        result.push_back(f->values());
    return result;
}
void exercise(const std::string& model, bool enabled, double temperature, int outers) {
    auto f=geometry();auto c=configuration(model,enabled,temperature);
    PS::EulerianEulerianPhaseSystem phases(c);phases.initialize(f);
    std::vector<PS::PhaseVectorField> previous;
    for(const auto& p:phases.phases()) previous.push_back(p.primitive.velocity);
    PS::PhaseSourceRegistry sources;
    if(model=="RPI") {
        WallHeatSetting wall;wall.patch="bottom";wall.heatFlux=7000;wall.wallTemperature=390;wall.coupleEnergy=true;
        sources.add(PS::makePhaseWallHeatSource({wall}));
    }
    FDM::SolverConfig solver;solver.numerics.timeRecipe=FDM::builtInTimeRecipe(FDM::TimeRecipeId::ForwardEuler);solver.numerics.maxDeltaT=1e-4;
    solver.pressure.coupling.preset=FDM::PressureCouplingPreset::PIMPLE;
    solver.pressure.coupling.outerCorrectors=outers;solver.pressure.coupling.pressureCorrectors=2;
    solver.pressure.coupling.nonOrthogonalCorrectors=1;
    System::BuildRequest request;request.templateOrigin=System::PhysicsTemplateKind::EulerianEulerian;
    request.phaseNames=c.eulerianEulerian.phaseNames;request.referencePhase="liquid";
    request.coupling=System::couplingRequestFrom(solver.pressure.coupling,true);
    const auto compiled=System::build(solver,request);
    Run::OpRegistry operations;std::map<std::string,int> count;
    Snapshot initial;
    for(const auto& p:phases.phases()) {
        initial.push_back(p.primary.phaseMass.values());initial.push_back(p.primary.phaseEnthalpy.values());
        for(int d=0;d<3;++d) initial.push_back(p.primary.momentum[d].values());
    }
    Snapshot baseline;double time=.3;int step=17;int probes=0;std::string preceding;
    for(const auto& id:compiled.runtime.report.requiredOperations) {
        operations.bind(id,"flow.eulerian-pressure",[&,id] {
            ++count[id];
            if(id==System::OpIds::EeDtCompute) {
                // stableTimeStep() evaluates the same sources at probe dt; this is
                // not an extra physical update and the following outer clears them.
                phases.computeInterphase(1e-4,previous);sources.assemble(phases,1e-4);++probes;
            }
            if(id==System::OpIds::EeInterphaseCompute) phases.computeInterphase(1e-4,previous);
            if(id==System::OpIds::EeSourcesAssemble) {
                require(preceding==System::OpIds::EeInterphaseCompute,"External source did not follow a fresh interphase assembly.");
                sources.assemble(phases,1e-4);
                const auto actual=sourceSnapshot(phases.sources());
                if(baseline.empty())baseline=actual;
                else require(actual==baseline,"Outer iteration accumulated source rates or wall diagnostics.");
                phases.sources().transferLedger.validate();
            }
            if(id==System::OpIds::EeTimeCommit) {time+=1e-4;++step;}
            preceding=id;
        });
    }
    // Only the real source operators execute here. Other numerical leaves are trace probes,
    // so this checks compiled placement and source budgets, not a complete CFD solution.
    for(int physicalStep=0;physicalStep<2;++physicalStep)
        Run::PlanExecutor::execute(compiled.solvePlan,operations);
    require(count[System::OpIds::EeInterphaseCompute]==outers*2 && count[System::OpIds::EeSourcesAssemble]==outers*2,
        "Source rate was assembled per phase/pressure/nonorth pass instead of per outer.");
    require(count[System::OpIds::EePressureSolve]==outers*8 && count[System::OpIds::EeStepCommit]==2
        && count[System::OpIds::EeTimeCommit]==2 && step==19 && probes==2,"Nested corrections repeated a physical commit.");
    near(time,.3002,"Source operation changed the physical clock.");
    Snapshot final;
    for(const auto& p:phases.phases()) {
        final.push_back(p.primary.phaseMass.values());final.push_back(p.primary.phaseEnthalpy.values());
        for(int d=0;d<3;++d) final.push_back(p.primary.momentum[d].values());
    }
    require(initial==final,"Source assembly advanced physical state.");
    const auto& s=phases.sources();const int cell=f.getIdx(5,4,3);
    double expectedMass=0,expectedHeat=0;
    if(model=="Lee" && enabled) {
        const auto& liquid=phases.phases()[0];const auto& gas=phases.phases()[1];
        expectedMass=Physics::PhaseChange::Lee::cellMassSource(.8,.2,temperature,373,
            liquid.primitive.density.values()[cell],gas.primitive.density.values()[cell],.01,.02);
    }
    if(model=="RPI") { expectedMass=enabled?.004:0;expectedHeat=14000; }
    near(s.mass[1].values()[cell],expectedMass,"Disabled RPI contributed mass, or Lee/RPI rate changed.");
    near(s.mass[0].values()[cell]+s.mass[1].values()[cell],0,"Mass exchange did not cancel.");
    for(int d=0;d<3;++d) near(s.momentum[0][d].values()[cell]+s.momentum[1][d].values()[cell],0,"Momentum exchange did not cancel.");
    near(s.energy[0].values()[cell]+s.energy[1].values()[cell],expectedHeat,"External heat or paired enthalpy counted twice.");
    if(model=="RPI") {
        near(s.wallBoilingMass.values()[cell],enabled?.004:0,"Inactive RPI has boiling diagnostics.");
        near(s.transferLedger.wallEnergyInput(cell),enabled?14000:0,"Wall ledger disagrees with active RPI contribution.");
        // Two inactive z planes share one physical thickness. Use one plane only.
        double power=0;
        for(int j=4;j<9;++j)for(int i=4;i<9;++i) {
            const int n=f.getIdx(i,j,3);
            power+=(s.energy[0].values()[n]+s.energy[1].values()[n])/f.Jac(i,j,3);
        }
        near(power,875,"Wall power differs from 7000 W/m2 times the supported 0.125 m2 patch measure.");
    }
    std::cout<<std::setprecision(17)<<model<<" enabled="<<enabled<<" T="<<temperature<<" outer="<<outers
        <<" mdot="<<s.mass[1].values()[cell]<<" qVol="<<expectedHeat<<" commits=2 time="<<time<<'\n';
}
}
int main() {
    for(int outer:{1,3}) {
        exercise("Lee",true,380,outer);exercise("Lee",true,360,outer);
        exercise("RPI",true,370,outer);exercise("RPI",false,370,outer);
    }
    auto bad=configuration("RPI",true,370);bad.phaseChange.selections.clear();
    bool rejected=false;
    try {PS::EulerianEulerianPhaseSystem invalid(bad);} catch(const std::runtime_error&) {rejected=true;}
    require(rejected,"Active RPI silently accepts missing closure selection.");
}
