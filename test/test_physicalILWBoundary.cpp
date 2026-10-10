#include "solver/boundary/reconstruction/SF_reconstruction.h"
#include "solver/boundary/reconstruction/ILW/SF_boundaryClosure.h"
#include "solver/boundary/SF_boundary.h"
#include "solver/boundary/SF_applicator.h"
#include "core/mesh/SF_dimension.h"
#include "methods/numerics/scalar/SF_scalarTransport.h"
#include "models/physics/fluidStateModel/SF_singleFluidStateModel.h"
#include "models/physics/EOS/SF_perfectGasEOS.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
using namespace SF;
namespace {
int checks=0;
void equal(double a,double b,const std::string& context,double tol=1.e-8) {
    ++checks;
    if(!std::isfinite(a) || std::abs(a-b)>tol*std::max(1.,std::abs(b)))
        throw std::runtime_error(context+": actual="+std::to_string(a)+" expected="+std::to_string(b));
}
void initialize(Field& f,bool threeD=false) {
    f.setup(15,15,threeD?15:2,3,6);
    for(int k=0;k<f.MZ();++k)for(int j=0;j<f.MY();++j)for(int i=0;i<f.MX();++i) {
        f.X(i,j,k)=.1*(i-3);f.Y(i,j,k)=.1*(j-3);f.Z(i,j,k)=.1*(k-3);f.Jac(i,j,k)=1000.;
        f(i,j,k,RHO)=1.2+.001*(i+j+k);
        f(i,j,k,RU)=2.*f(i,j,k,RHO);f(i,j,k,RV)=.3*f(i,j,k,RHO);f(i,j,k,RW)=.4*f(i,j,k,RHO);
        Boundary::setEnergyFromPressure(f,i,j,k,101500.);
        f(i,j,k,5)=7.;
    }
}
void markStale(Field& f,int axis,int side,int var,double value) {
    const int at=side<0?3:17;
    for(int k=3;k<3+f.NZ();++k)for(int j=3;j<18;++j)for(int i=3;i<18;++i) {
        const int p[3]={i,j,k};if(p[axis]==at)f(i,j,k,var)=value;
    }
}
void basic(int order,int axis,int side,bool threeD=false) {
    Field f;initialize(f,threeD);
    int p[3]={10,10,threeD?10:3};p[axis]=side<0?3:17;
    auto recipe=Boundary::Reconstruction::Selection{Boundary::Reconstruction::Kind::ILW,order};
    markStale(f,axis,side,5,-99.);
    Boundary::Reconstruction::applyScalar(f,p[0],p[1],p[2],axis,5,0.,Boundary::Law::Kind::ZeroGradient,recipe);
    equal(f(p[0],p[1],p[2],5),7.,"Neumann real point must be refreshed");
    for(int layer=1;layer<=3;++layer) {
        int g[3]={p[0],p[1],p[2]};g[axis]+=side*layer;
        equal(f(g[0],g[1],g[2],5),7.,"Neumann ghost");
        f(g[0],g[1],g[2],5)=-123.;
        Boundary::Reconstruction::applyScalar(f,g[0],g[1],g[2],axis,5,0.,Boundary::Law::Kind::ZeroGradient,recipe);
        equal(f(g[0],g[1],g[2],5),7.,"direct ghost entry");
    }
    // U has a Neumann condition; a varying density must not turn it into one on rhoU.
    markStale(f,axis,side,RU,0.);markStale(f,axis,side,RV,0.);markStale(f,axis,side,RW,0.);
    Boundary::Reconstruction::applyVector(f,p[0],p[1],p[2],axis,RU,{},Boundary::Law::Kind::ZeroGradient,recipe);
    for(int comp=0;comp<3;++comp)equal(f(p[0],p[1],p[2],RU+comp)/f(p[0],p[1],p[2],RHO),comp==0?2.:comp==1?.3:.4,"primitive velocity Neumann");
    Boundary::Reconstruction::applyVector(f,p[0],p[1],p[2],axis,RU,{1.,.6,.2},Boundary::Law::Kind::FixedValue,recipe);
    const double fixed[3]={1.,.6,.2};
    for(int comp=0;comp<3;++comp)equal(f(p[0],p[1],p[2],RU+comp)/f(p[0],p[1],p[2],RHO),fixed[comp],"fixed velocity");
    Boundary::Reconstruction::applyVector(f,p[0],p[1],p[2],axis,RU,{},Boundary::Law::Kind::Symmetry,recipe);
    for(int comp=0;comp<3;++comp)equal(f(p[0],p[1],p[2],RU+comp)/f(p[0],p[1],p[2],RHO),comp==axis?0.:comp==0?2.:comp==1?.3:.4,"symmetry normal and current tangents");
    // Recover pressure from current interior after poisoning the whole real face.
    const int at=p[axis];
    for(int k=3;k<3+f.NZ();++k)for(int j=3;j<18;++j)for(int i=3;i<18;++i) {
        const int q[3]={i,j,k};if(q[axis]==at)Boundary::setEnergyFromPressure(f,i,j,k,101325.);
    }
    const std::string name=axis==0?"Right":axis==1?"Top":"Back";
    f.setBoundarySets({{name,{f.getIdx(p[0],p[1],p[2])}}});
    std::vector<BCSetting<double>> bc{{name,ZERO_GRADIENT,0.}};
    Boundary::updateEnergyFromPressure(f,bc,true,order);
    equal(Boundary::pressureAt(f,p[0],p[1],p[2]),101500.,"pressure Neumann real");
    bc[0].type=FIXED_VALUE;bc[0].value=101000.;
    Boundary::updateEnergyFromPressure(f,bc,true,order);
    equal(Boundary::pressureAt(f,p[0],p[1],p[2]),101000.,"fixed pressure");
}
void thermal(int order) {
    Field f;initialize(f);
    for(int k=0;k<f.MZ();++k)for(int j=0;j<f.MY();++j)for(int i=0;i<f.MX();++i)
        Boundary::setEnergyFromTemperature(f,i,j,k,300.);
    f.setBoundarySets({{"Top",{f.getIdx(10,17,3)}}});
    for(auto type:{ThermalBCType::ZeroGradient,ThermalBCType::Adiabatic,ThermalBCType::FixedTemperature,ThermalBCType::HeatFlux}) {
        Boundary::setEnergyFromTemperature(f,10,17,3,280.);
        const double kappa=.001*DefaultIdealGasGamma*DefaultIdealGasConstant/(DefaultIdealGasGamma-1.)/.72;
        std::vector<ThermalBCSetting> bc{{"Top",type,type==ThermalBCType::HeatFlux?2.*kappa:310.}};
        Boundary::updateEnergyFromThermalBoundary(f,bc,.001,.72,true,order);
        equal(Boundary::temperatureAt(f,10,17,3),type==ThermalBCType::FixedTemperature?310.:type==ThermalBCType::HeatFlux?300.2:300.,"thermal real");
        for(int layer=1;layer<=3;++layer)
            equal(Boundary::temperatureAt(f,10,17+layer,3),type==ThermalBCType::FixedTemperature?310.+10.*layer:type==ThermalBCType::HeatFlux?300.2+.2*layer:300.,"thermal ghost gradient/value");
    }
}
void cornerAndParity(int order) {
    Field f;initialize(f);
    auto recipe=Boundary::Reconstruction::Selection{Boundary::Reconstruction::Kind::ILW,order};
    markStale(f,1,1,5,-99.);
    Boundary::Reconstruction::applyScalar(f,3,17,3,1,5,0.,Boundary::Law::Kind::ZeroGradient,recipe);
    equal(f(3,17,3,5),7.,"corner uses explicit y normal");
    equal(f(3,18,3,5),7.,"corner y ghost");
    for(int k=0;k<f.MZ();++k)for(int j=0;j<f.MY();++j)for(int i=0;i<f.MX();++i) {
        const double s=.1*(j-17);
        f(i,j,k,RV)=s*f(i,j,k,RHO);
    }
    Boundary::Reconstruction::applyVector(f,10,17,3,1,RU,{},Boundary::Law::Kind::Symmetry,recipe);
    equal(f(10,17,3,RV),0.,"symmetry zero normal");
    if(f(10,18,3,RV)<=0.)throw std::runtime_error("symmetry ghost must preserve odd normal extension");
    // Higher derivatives must remain active, not be replaced by constant copying.
    for(int k=0;k<f.MZ();++k)for(int j=0;j<f.MY();++j)for(int i=0;i<f.MX();++i)f(i,j,k,5)=7.+.01*std::pow(.1*(j-17),2);
    Boundary::Reconstruction::applyScalar(f,10,17,3,1,5,0.,Boundary::Law::Kind::ZeroGradient,recipe);
    if(std::abs(f(10,20,3,5)-f(10,17,3,5))<1.e-7)throw std::runtime_error("Neumann was silently reduced to constant copying");
}
void scalarStorage(int order) {
    Field f;initialize(f);ScalarField c;c.setupLike(f,"C",7.);
    for(int i=3;i<18;++i)c(i,17,3)=-99.;
    f.setBoundarySets({{"Top",{f.getIdx(10,17,3)}}});
    for(auto law:{ZERO_GRADIENT,SYMMETRY,FIXED_VALUE}) {
        FDM::Scalar::applyBoundaryConditions(f,c,{{"Top",law,8.}},true,order);
        equal(c(10,17,3),law==FIXED_VALUE?8.:7.,"independent scalar real");
        equal(c(10,18,3),law==FIXED_VALUE?9.:7.,"independent scalar ghost");
    }
}
void boundThermodynamics() {
    Field f;f.setup(15,15,2,3,5);
    auto eos=std::make_shared<Physics::EOS::PerfectGasEOS>(1.67,208.1,1234.);
    f.setStateModel(std::make_shared<Physics::FluidStateModel::SingleFluidStateModel>(eos,.001,2.));
    for(int k=0;k<f.MZ();++k)for(int j=0;j<f.MY();++j)for(int i=0;i<f.MX();++i) {
        f.X(i,j,k)=.1*i;f.Y(i,j,k)=.1*j;f.Z(i,j,k)=.1*k;
        f(i,j,k,RHO)=1.2;f(i,j,k,RU)=2.4;f(i,j,k,RV)=.36;f(i,j,k,RW)=0.;
        Boundary::setEnergyFromTemperature(f,i,j,k,300.);
    }
    f.setBoundarySets({{"Top",{f.getIdx(10,17,3)}}});
    Boundary::updateEnergyFromThermalBoundary(f,{{"Top",ThermalBCType::HeatFlux,4.}},0.,.72,true,5);
    f.invalidateThermodynamicCache(); // Match the existing conservative-state publication contract.
    equal(Boundary::temperatureAt(f,10,17,3),300.2,"active EOS thermal real");
    equal(Boundary::temperatureAt(f,10,18,3),300.4,"active EOS conductivity");
    Boundary::updateEnergyFromPressure(f,{{"Top",FIXED_VALUE,101000.}},true,5);
    f.invalidateThermodynamicCache();
    equal(Boundary::pressureAt(f,10,17,3),101000.,"active EOS pressure inversion");
    const double ke=.5*(2.4*2.4+.36*.36)/1.2;
    equal(f(10,17,3,E),101000./.67+1.2*1234.+ke,"active EOS reference energy and kinetic energy");
}
void linearSymmetry() {
    for(int axis:{0,1})for(int side:{-1,1}) {
        Field f;initialize(f);
        for(int k=0;k<f.MZ();++k)for(int j=0;j<f.MY();++j)for(int i=0;i<f.MX();++i) {
            f(i,j,k,RHO)=1.;f(i,j,k,RU)=2.;f(i,j,k,RV)=.3;f(i,j,k,RW)=.4;
        }
        for(int c=0;c<3;++c)markStale(f,axis,side,RU+c,0.);
        int p[3]={10,10,3};p[axis]=side<0?3:17;
        Boundary::Reconstruction::applyVector(f,p[0],p[1],p[2],axis,RU,{},Boundary::Law::Kind::Symmetry,
            {Boundary::Reconstruction::Kind::Linear,0});
        for(int c=0;c<3;++c)equal(f(p[0],p[1],p[2],RU+c),c==axis?0.:c==0?2.:c==1?.3:.4,"linear symmetry current tangents");
    }
}
void orderedThermalPublication(int order) {
    Field f;f.setup(15,15,2,3,5);
    auto eos=std::make_shared<Physics::EOS::PerfectGasEOS>(1.4,287.05);
    f.setStateModel(std::make_shared<Physics::FluidStateModel::SingleFluidStateModel>(eos,.001,2.));
    for(int k=0;k<f.MZ();++k)for(int j=0;j<f.MY();++j)for(int i=0;i<f.MX();++i) {
        f.X(i,j,k)=.1*i;f.Y(i,j,k)=.1*j;f.Z(i,j,k)=.1*k;
        f(i,j,k,RHO)=1.2;f(i,j,k,RU)=2.4;f(i,j,k,RV)=0.;f(i,j,k,RW)=0.;
        Boundary::setEnergyFromTemperature(f,i,j,k,300.);
    }
    f.setBoundarySets({{"Top",{f.getIdx(16,17,3)}},
        {"Right",{f.getIdx(17,16,3),f.getIdx(17,17,3)}}});
    equal(Boundary::temperatureAt(f,16,17,3),300.,"prime old EOS cache at an intersecting face");
    Boundary::updateEnergyFromThermalBoundary(f,{{"Top",ThermalBCType::FixedTemperature,310.},
        {"Right",ThermalBCType::ZeroGradient,0.}},.001,.72,order!=0,order);
    equal(Boundary::temperatureAt(f,16,17,3),310.,"thermal writer invalidates its old EOS view");
    equal(Boundary::temperatureAt(f,17,17,3),310.,"next thermal patch consumes current published energy");
    equal(Boundary::temperatureAt(f,18,17,3),310.,"corner ghost uses current thermal donor");
}
void empty() {
    Field f;initialize(f);
    Boundary::Reconstruction::Selection recipe{Boundary::Reconstruction::Kind::ILW,9};
    for(int k=0;k<3;++k)f(10,10,k,5)=-99.;
    Boundary::Reconstruction::applyScalar(f,10,10,3,2,5,0.,Boundary::Law::Kind::Empty,recipe);
    for(int k=0;k<3;++k)equal(f(10,10,k,5),7.,"empty inactive-axis copy");
}
void wallLaws(int order,int axis,int side,bool threeD=false) {
    using namespace Boundary;
    auto recipe=Reconstruction::fromILWSetting(order!=0,order);
    for(auto law:{Law::Kind::Slip,Law::Kind::NoSlip}) {
        Field f;initialize(f,threeD);
        int p[3]={10,10,threeD?10:3};p[axis]=side<0?3:17;
        for(int c=0;c<3;++c)markStale(f,axis,side,RU+c,0.);
        Reconstruction::applyVector(f,p[0],p[1],p[2],axis,RU,{},law,recipe);
        for(int layer=0;layer<=3;++layer) {
            int g[3]={p[0],p[1],p[2]};g[axis]+=side*layer;
            for(int c=0;c<3;++c) {
                const double interior=c==0?2.:c==1?.3:.4;
                const double expected=law==Law::Kind::NoSlip ? 0. : c==axis
                    ? (layer>0 && order!=0 ? -interior*layer : 0.) : interior;
                if(law==Law::Kind::NoSlip && layer>0 && order!=0) continue;
                equal(f(g[0],g[1],g[2],RU+c)/f(g[0],g[1],g[2],RHO),expected,"wall primitive velocity and independent density");
            }
        }
        Field neutral;initialize(neutral,threeD);
        for(int c=0;c<3;++c)markStale(neutral,axis,side,RU+c,0.);
        Reconstruction::applyVector(neutral,p[0],p[1],p[2],axis,RU,{},
            law==Law::Kind::NoSlip?Law::Kind::FixedValue:Law::Kind::Symmetry,recipe);
        for(int layer=0;layer<=3;++layer) {
            int g[3]={p[0],p[1],p[2]};g[axis]+=side*layer;
            for(int c=0;c<3;++c)equal(f(g[0],g[1],g[2],RU+c),neutral(g[0],g[1],g[2],RU+c),"wall and expanded neutral constraint are identical");
        }
        // Velocity choice cannot replace the separately prescribed pressure or temperature.
        const std::string name=axis==0?"Right":axis==1?"Top":"Back";
        f.setBoundarySets({{name,{f.getIdx(p[0],p[1],p[2])}}});
        updateEnergyFromPressure(f,{{name,FIXED_VALUE,101000.}},order!=0,order);
        equal(pressureAt(f,p[0],p[1],p[2]),101000.,"wall does not force pressure Neumann");
        updateEnergyFromThermalBoundary(f,{{name,ThermalBCType::FixedTemperature,305.}},.001,.72,order!=0,order);
        equal(temperatureAt(f,p[0],p[1],p[2]),305.,"wall does not force adiabatic temperature");
    }
    Field f;initialize(f,threeD);int p[3]={10,10,threeD?10:3};p[axis]=side<0?3:17;
    for(auto law:{Law::Kind::Slip,Law::Kind::NoSlip}) {
        bool rejected=false;try {Reconstruction::applyScalar(f,p[0],p[1],p[2],axis,5,0.,law,recipe);}catch(const std::runtime_error&){rejected=true;}
        equal(rejected,1.,"wall law rejects scalar target");
        rejected=false;try {Reconstruction::applyVector(f,p[0],p[1],p[2],axis,RU,{1.,0.,0.},law,recipe);}catch(const std::runtime_error&){rejected=true;}
        equal(rejected,1.,"stationary wall rejects conflicting moving velocity");
    }
    f.setBoundarySets({{"RemotePhysicalWall",{}}});
    const double before=f(10,10,threeD?10:3,E);
    updateEnergyFromThermalBoundary(f,{{"RemotePhysicalWall",ThermalBCType::FixedTemperature,305.}},.001,.72,order!=0,order);
    equal(f(10,10,threeD?10:3,E),before,"thermal boundary with zero local partition members has no write");
}
}
int main(int argc,char**) {
    if(argc>1) {
        Math::resetActiveDirections();Math::deactivateDirection(2);
        Field f;f.setup(4,4,2,3,5);
        for(int k=0;k<f.MZ();++k)for(int j=0;j<f.MY();++j)for(int i=0;i<f.MX();++i) {f.X(i,j,k)=i;f.Y(i,j,k)=j;f.Z(i,j,k)=k;f(i,j,k,0)=1.;}
        Boundary::Reconstruction::applyScalar(f,3,4,3,0,0,0.,Boundary::Law::Kind::ZeroGradient,{Boundary::Reconstruction::Kind::ILW,9});
        return 2;
    }
    try {
        Math::resetActiveDirections();Math::deactivateDirection(2);
        for(int order:{3,5,7,9}) {
            for(int axis:{0,1})for(int side:{-1,1})basic(order,axis,side);
            thermal(order);cornerAndParity(order);scalarStorage(order);
            for(int axis:{0,1})for(int side:{-1,1})wallLaws(order,axis,side);
            std::cout<<"ILW"<<order<<" scalar/vector/pressure/thermal/corner checks passed\n";
        }
        boundThermodynamics();linearSymmetry();empty();Math::resetActiveDirections();
        for(int order:{3,5,7,9}) {basic(order,2,order<7?-1:1,true);wallLaws(order,2,order<7?-1:1,true);}
        Math::deactivateDirection(2);
        for(int axis:{0,1})for(int side:{-1,1})wallLaws(0,axis,side);
        for(int order:{0,3,5,7,9})orderedThermalPublication(order);
        std::cout<<checks<<" physical ILW assertions passed\n";
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
