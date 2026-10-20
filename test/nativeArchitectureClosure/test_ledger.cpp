#include "models/physics/phaseChange/Lee/SF_Lee.h"
#include "models/physics/phaseChange/RPI/SF_heatFluxPartition.h"
#include "models/physics/interphaseTransfer/SF_transferLedger.h"
#include "models/physics/fluidStateModel/SF_homogeneousMultiphaseStateModel.h"
#include <iostream>
#include <stdexcept>
using namespace SF;
void require(bool v,const char* s){if(!v)throw std::runtime_error(s);}
int main(){
 namespace PC=Physics::PhaseChange;namespace MP=Physics::Multiphase;namespace IT=Physics::InterphaseTransfer;
 Field f;f.setup(1,1,1,1,5);ScalarField alpha,t;alpha.setupLike(f,"alpha");t.setupLike(f,"T");
 MP::MultiPhaseConfig c;c.phaseChange.latentHeat=2e6;
 c.phaseChange.selections={{"convectiveheatfluxmodel","fixed"},{"quenchingheatfluxmodel","fixed"}};
 c.phaseChange.coefficients={{"convectiveheatflux",1000},{"quenchingheatflux",2000}};
 MP::PhaseProperties liquid,gas;liquid.name="liquid";gas.name="gas";
 PC::ModelContext context{f,alpha,t,c,liquid,gas,liquid};
 double rate=PC::Lee::cellMassSource(.8,.2,380,373,1000,1,2,3);
 require(rate==2*.8*1000*(380.-373)/373,"Lee evaporation law changed");
 double condensation=PC::Lee::cellMassSource(.8,.2,360,373,1000,1,2,3);
 require(condensation==-3*.2*1*(373.-360)/373,"Lee condensation sign changed");
 const auto heat=PC::RPI::evaluateHeatFluxPartition(context,390,370,.002);
 require(heat.convective==1000 && heat.quenching==2000 && heat.evaporative==4000 && heat.total==7000,"RPI wall budget changed");
 IT::TransferLedger ledger(1,2);
 ledger.addInternalMassTransfer(0,0,1,rate);
 ledger.addInternalMomentumTransfer(0,0,1,{2*rate,0,0});
 ledger.addInternalEnergyTransfer(0,0,1,rate*1e5);
 ledger.addWallEnergyInput(0,0,heat.convective+heat.quenching);
 ledger.addWallEnergyInput(0,1,heat.evaporative);

 ledger.validate();
 require(ledger.phaseMassSource(0,0)+ledger.phaseMassSource(0,1)==0,"Internal mass imbalance");
 require(ledger.phaseMomentumSource(0,0,0)+ledger.phaseMomentumSource(0,1,0)==0,"Internal momentum imbalance");
 require(ledger.phaseEnergySource(0,0)+ledger.phaseEnergySource(0,1)==7000,"Wall/work energy counted incorrectly");
 require(ledger.mixtureEnergySource(0)==7000 && ledger.mechanicalEnergyConversion(0)==0,"External energy or mechanical conversion counted twice");
 using Physics::FluidStateModel::Component;using Physics::FluidStateModel::HomogeneousMultiphaseStateModel;
 auto eos=std::make_shared<Physics::EOS::StiffenedGasEOS>(1.4,718,0,0);
 HomogeneousMultiphaseStateModel model({Component{"liquid","bulk",eos},Component{"gas","bulk",eos}});
 IT::SourceBundle sources(1,model.variableCount());model.mapTransfersToSources(ledger,{0,1},sources);
 require(sources(0,0)+sources(0,1)==0 && sources(0,model.energyIndex())==7000,"Conservative energy contains duplicate latent or internal exchange");
 auto vapor=std::make_shared<Physics::EOS::StiffenedGasEOS>(1.4,718,0,2e6);
 HomogeneousMultiphaseStateModel latent({Component{"liquid","bulk",eos},Component{"gas","bulk",vapor}});
 auto q=latent.conservativeFromPressureTemperature({.8,.2},101325,300,{2,0,0});
 const auto initial=latent.close(q.data(),latent.variableCount());const double initialEnergy=q[latent.energyIndex()];
 IT::TransferLedger internal(1,2);internal.addInternalMassTransfer(0,0,1,.001);
 IT::SourceBundle source(1,latent.variableCount());latent.mapTransfersToSources(internal,{0,1},source);
 for(int n=0;n<latent.variableCount();++n)q[n]+=source(0,n)*.01;
 const auto cooled=latent.close(q.data(),latent.variableCount());
 require(q[latent.energyIndex()]==initialEnergy && std::abs(cooled.density-initial.density)<1e-12
     && cooled.temperature<initial.temperature,"EOS reference latent energy added a second total-energy source");
 ledger.addInternalMechanicalEnergy(0,0,30);ledger.validate();
 require(ledger.phaseEnergySource(0,0)+ledger.phaseEnergySource(0,1)==7030,"Internal mechanical energy partition lost");
 bool rejected=false;try{IT::SourceBundle bad(1,model.variableCount());model.mapTransfersToSources(ledger,{0,1},bad);}catch(const std::runtime_error&){rejected=true;}
 require(rejected,"Homogeneous silently accepts unrepresented phase mechanical conversion");
 std::cout<<"Lee signed rate, RPI 7000 W/m2 partition, paired mass/momentum/enthalpy, external wall/work ledger, homogeneous conservative mapping passed; no CFD execution claim\n";
}
