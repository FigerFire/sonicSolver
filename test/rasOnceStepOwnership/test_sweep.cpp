#include "solver/algorithm/transport/SF_transportBinding.h"
#include "models/turbulence/SF_turbulence.h"
#include "models/physics/fluidStateModel/SF_factory.h"
#include "core/mesh/SF_dimension.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
using namespace SF;
void dump(const char* phase,int step,const Turbulence::Manager& m) {
 std::cout<<phase<<' '<<step<<'\n'<<std::hexfloat;
 for(auto slot:{Turbulence::ScalarSlot::K,Turbulence::ScalarSlot::Epsilon,Turbulence::ScalarSlot::Omega,Turbulence::ScalarSlot::EddyMu}) {
  for(double x:m.scalarFields().values(slot))std::cout<<x<<' ';std::cout<<'\n';
 }
}
int main(){
 for(auto kind:{FDM::TurbulenceModelKind::kEpsilon,FDM::TurbulenceModelKind::kOmegaSST,FDM::TurbulenceModelKind::Smagorinsky}) {
  Field f;f.setup(9,5,2,3,5);f.setStateModel(Physics::FluidStateModel::makeSingleFluidPerfectGas(1.4,287.05,.001,.72));
  Math::resetActiveDirections();Math::deactivateDirection(2);
  FDM::TurbulenceConfig c;c.enabled=true;c.model=kind;c.family=kind==FDM::TurbulenceModelKind::Smagorinsky?FDM::TurbulenceFamily::LES:FDM::TurbulenceFamily::RAS;c.laminarDynamicViscosity=.001;
  for(int n=0;n<f.TotalSize();++n){int i,j,k;f.getIJK(n,i,j,k);double x=(i-3)*.2,y=(j-3)*.3;
   f.X(i,j,k)=x;f.Y(i,j,k)=y;f.Z(i,j,k)=(k-3)*.1;f.Jac(i,j,k)=1/(.2*.3*.1);f.XiX(i,j,k)=5;f.EtY(i,j,k)=1/.3;f.ZeZ(i,j,k)=10;f.wallDistance(i,j,k)=.08+.01*(i+1);
   double rho=1.2+.02*std::sin(x),u=.3+.04*x*y,v=.02*x;
   f(i,j,k,0)=rho;f(i,j,k,1)=rho*u;f(i,j,k,2)=rho*v;f(i,j,k,3)=0;f(i,j,k,4)=100000/.4+.5*rho*(u*u+v*v);f.CellFlag(i,j,k)=FLUID_CELL;
   std::string set="point"+std::to_string(n);f.getAllSets()[set]={n};
   c.scalars.kInitial.push_back({set,FIXED_VALUE,.02+.001*x*x+.0005*y});
   c.scalars.epsilonInitial.push_back({set,FIXED_VALUE,.03+.002*x+.001*y*y});
   c.scalars.omegaInitial.push_back({set,FIXED_VALUE,2+.05*x+.03*y*y});
   if(i>=3&&i<=11&&j>=3&&j<=7&&k>=3&&k<=4){
    if(i==3)f.getAllSets()["left"].push_back(n);if(i==11)f.getAllSets()["right"].push_back(n);
    if(j==3)f.getAllSets()["bottom"].push_back(n);if(j==7)f.getAllSets()["top"].push_back(n);
    if(k==3||k==4)f.getAllSets()["empty"].push_back(n);
   }
  }
  c.scalars.kBoundary={{"left",FIXED_VALUE,.015},{"right",ZERO_GRADIENT,0},{"bottom",SYMMETRY,0},{"top",FIXED_VALUE,.025},{"empty",EMPTY,0}};
  c.scalars.epsilonBoundary=c.scalars.kBoundary;c.scalars.epsilonBoundary[0].value=.03;c.scalars.epsilonBoundary[3].value=.04;
  c.scalars.omegaBoundary=c.scalars.kBoundary;c.scalars.omegaBoundary[0].value=1.5;c.scalars.omegaBoundary[3].value=2.5;
  Turbulence::Manager manager(c);if(!manager.initialize(f))throw std::runtime_error("Missing model");
  SolverAlgorithm::initializeTransport(f,manager);
  if(c.family==FDM::TurbulenceFamily::RAS) {
   const auto port=manager.onceStepTransport();
   const auto slot=kind==FDM::TurbulenceModelKind::kEpsilon?Turbulence::ScalarSlot::Epsilon:Turbulence::ScalarSlot::Omega;
   if(port.variables[0]!=&manager.scalarFields().values(Turbulence::ScalarSlot::K)
     ||port.variables[1]!=&manager.scalarFields().values(slot)
     ||port.eddyViscosity!=&manager.scalarFields().values(Turbulence::ScalarSlot::EddyMu))
    throw std::runtime_error("Once-step port replaced the original physical arrays");
  }
  std::cout<<"model "<<int(kind)<<'\n';dump("initial",0,manager);
  for(int step=0;step<5;++step){SolverAlgorithm::applyTransportBoundary(f,manager);dump("boundary-ready",step,manager);SolverAlgorithm::updateOnceStepTransport(f,manager,.00002);dump("updated",step,manager);SolverAlgorithm::applyTransportBoundary(f,manager);dump("published",step,manager);}
 }
}
