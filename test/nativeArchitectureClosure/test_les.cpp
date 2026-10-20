#include "models/turbulence/LES/SF_LES.h"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <iostream>
using namespace SF;
int main(){
 Field f;f.setup(5,5,1,1,5);Turbulence::ScalarFields state;
 FDM::TurbulenceConfig config;Turbulence::LES::SmagorinskyModel model;
 for(int k=0;k<f.MZ();++k)for(int j=0;j<f.MY();++j)for(int i=0;i<f.MX();++i){
  f.X(i,j,k)=i;f.Y(i,j,k)=j;f.Z(i,j,k)=k;
  f.XiX(i,j,k)=f.EtY(i,j,k)=f.ZeZ(i,j,k)=f.Jac(i,j,k)=1;
  f(i,j,k,0)=1;f(i,j,k,1)=j;f(i,j,k,4)=100+j*j*.5;
 }
 const int i=3,j=3,k=1;
 for(double rho:{1.,1000.,1e-15}){
  for(int cell=0;cell<f.TotalSize();++cell){int a,b,c;f.getIJK(cell,a,b,c);f(a,b,c,0)=rho;f(a,b,c,1)=rho*b;f(a,b,c,4)=rho*(100+b*b*.5);}
  double d=Turbulence::characteristicFilterWidth(f,i,j,k,config.coefficients.filterScale);
  double s=Turbulence::strainRateMagnitude(f,i,j,k),cs=config.coefficients.cSmagorinsky;
  double expected=rho*(cs*d)*(cs*d)*s,actual=model.eddyDynamicViscosity(f,state,config,i,j,k);
  if(actual!=expected || !std::isfinite(actual))throw std::runtime_error("Smagorinsky valid-rho formula changed");
 }
 for(double rho:{0.,-1.,std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()}){
  f(i,j,k,0)=rho;bool rejected=false;
  try{model.eddyDynamicViscosity(f,state,config,i,j,k);}catch(const std::runtime_error&){rejected=true;}
  if(!rejected)throw std::runtime_error("LES invalid density hidden clamp");
 }
 std::cout<<"valid density direct formula (including 1e-15); zero/negative/NaN/Inf rejected\n";
}
