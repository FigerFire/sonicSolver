#include "models/ibm/method/ghost/SF_viscousWall.h"
#include "models/ibm/SF_ibmSystemContribution.h"
#include "models/ibm/descriptor/SF_algorithmDescriptor.h"
#include "models/turbulence/SF_turbulenceSystemContribution.h"
#include "solver/system/SF_systemBuilder.h"
#include "solver/system/SF_stateRealizer.h"
#include "models/turbulence/SF_turbulence.h"
#include "methods/numerics/viscous/SF_viscous.h"
#include <iostream>
using namespace SF;
void require(bool ok,const char* message) {if(!ok)throw std::runtime_error(message);}
int main() {
    try {
        Math::resetActiveDirections();Math::deactivateDirection(2);
        double previousGhost=0,previousShear=0,previousCentral=0;
        for (double h:{.04,.02,.01,.005}) {
            Field field;field.setup(7,7,2,1,5);field.clearCellFlags(FLUID_CELL);
            IBM::IBMGeometry geometry;geometry.setup(field.MX(),field.MY(),field.MZ());
            IBM::GhostILW::Storage storage;storage.setup(field.TotalSize());
            std::vector<int> cells;
            for(int k=1;k<=2;++k)for(int j=1;j<=7;++j)for(int i=1;i<=7;++i) {
                const double x=(i-4)*h,y=(j-2)*h;
                field.X(i,j,k)=x;field.Y(i,j,k)=y;field.Z(i,j,k)=.1*(k-1);
                field.setWallDistance(i,j,k,std::abs(y));
                field.XiX(i,j,k)=field.EtY(i,j,k)=1/h;
                if (j<3)field.CellFlag(i,j,k)=IBM_GHOST_CELL;
                const double u=y+2*y*y+.3*y*x+.7*y*y*y,rho=1+.2*x,p=10+.1*x*x;
                field(i,j,k,RHO)=rho;field(i,j,k,RU)=rho*u;field(i,j,k,RV)=field(i,j,k,RW)=0;
                field(i,j,k,E)=p/.4+.5*rho*u*u;
                if (j>=3 && j<=5)cells.push_back(field.getIdx(i,j,k));
            }
            // Evaluate exactly one ghost at y=-h, wall y=0; no source geometry fallback.
            const int i=4,j=1,k=1,ghost=field.getIdx(i,j,k);
            geometry.signedDistance(i,j,k)=-h;
            geometry.setGeometry(i,j,k,{0,0,0},{0,h,0},{0,1,0},{0,0,0});
            storage.setFluidSamples(ghost,cells);
            const auto plan=IBM::GhostIBM::buildViscousWallPlan(field,storage,geometry,i,j,k);
            double derivative=0;
            for(std::size_t q=0;q<cells.size();++q) {
                int a,b,c;field.getIJK(cells[q],a,b,c);derivative+=plan.normalDerivative[q]*field(a,b,c,RU)/field(a,b,c,RHO);
            }
            IBM::GhostIBM::applyViscousWall(field,{plan},1.4);
            const double u=field(i,j,k,RU)/field(i,j,k,RHO);
            geometry.setGeometry(i,j,k,{-h,0,0},{0,h,0},{0,1,0},{0,0,0});
            const auto shifted=IBM::GhostIBM::buildViscousWallPlan(field,storage,geometry,i,j,k);
            IBM::GhostIBM::applyViscousWall(field,{shifted},1.4);
            require(std::abs(field(i,j,k,RU)/field(i,j,k,RHO)-u)<1e-12,
                "Moving a planar wall anchor tangentially changed the physical ghost evaluation.");
            const double ghostError=std::abs(u-(-h+2*h*h-.7*h*h*h)),shearError=std::abs(derivative-1);
            require(std::abs(field(i,j,k,RHO)-1)<1e-12,"Neumann density was not reproduced.");
            require(ghostError<7*h*h*h && shearError<4*h*h,"Viscous boundary failed manufactured wall accuracy.");
            if(previousGhost)require(previousGhost/ghostError>7.9 && previousShear/shearError>3.9,"Viscous wall/shear lost grid convergence.");
            std::cout<<"h="<<h<<" ghostError="<<ghostError<<" shearError="<<shearError<<'\n';
            previousGhost=ghostError;previousShear=shearError;
            // The production Central2/Newtonian kernel consumes the same reconstructed wall stencil.
            geometry.signedDistance(4,2,1)=0;
            geometry.setGeometry(4,2,1,{0,0,0},{0,0,0},{0,1,0},{0,0,0});
            storage.setFluidSamples(field.getIdx(4,2,1),cells);
            const auto surface=IBM::GhostIBM::buildViscousWallPlan(field,storage,geometry,4,2,1);
            IBM::GhostIBM::applyViscousWall(field,{surface},1.4);
            require(field(4,2,1,RU)==0,"On-surface node is not exactly no-slip.");
            Viscous::TransportProperties transport;transport.mu=.001;
            const auto gradient=Viscous::velocityGradientAt(field,4,2,1,Viscous::CentralOrder::SECOND,transport);
            const auto stress=Constitutive::Newtonian::stress(transport.mu,gradient);
            const double centralError=std::abs(gradient[0].y-1);
            require(centralError<5*h*h,"Production Central2 wall gradient differs from analytic shear.");
            require(std::abs(stress.xy-transport.mu)<transport.mu*5*h*h,"Production Newtonian wall stress differs from analytic shear.");
            if(previousCentral)require(previousCentral/centralError>3.9,"Production wall shear lost second-order convergence.");
            previousCentral=centralError;
            std::cout<<"Central2 shearError="<<centralError<<'\n';
            require(std::abs(Viscous::physicalGradientAt(field,4,2,1,3,
                Viscous::CentralOrder::SECOND,transport).y)<1e-10,"Adiabatic wall has nonzero normal thermal gradient.");
            std::vector<double> kValues(field.TotalSize()),omega(field.TotalSize()),muT(field.TotalSize(),.01);
            const double wallOmega=60*.00001/(.075*h*h);
            for (int cell:cells) {
                int a,b,c;field.getIJK(cell,a,b,c);const double y=field.Y(a,b,c);
                field(a,b,c,RHO)=1;
                kValues[cell]=.001*y*y;omega[cell]=wallOmega+y+y*y;
            }
            IBM::GhostIBM::applySSTWall(field,{plan},.00001,kValues,omega,muT);
            require(std::abs(kValues[ghost]-.001*h*h)<1e-12,"SST k wall profile differs from analytic Dirichlet reference.");
            require(std::abs(omega[ghost]-(wallOmega-h+h*h))<1e-10,"SST omega wall profile differs from analytic TMR boundary reference.");
            require(muT[ghost]==0,"SST immersed ghost retained eddy viscosity.");

        }
        struct Wall final:FDM::IImmersedTurbulenceBoundary {
            mutable int calls=0;
            void applySST(const Field&,double,std::vector<double>&,std::vector<double>&,std::vector<double>&) const override {++calls;}
        } boundary;
        Field frozen;frozen.setup(5,5,2,1,5);frozen.clearCellFlags(FLUID_CELL);
        for(int k=0;k<frozen.MZ();++k)for(int j=0;j<frozen.MY();++j)for(int i=0;i<frozen.MX();++i) {
            frozen(i,j,k,RHO)=1;frozen(i,j,k,RU)=frozen(i,j,k,RV)=frozen(i,j,k,RW)=0;frozen(i,j,k,E)=25;
            frozen.setWallDistance(i,j,k,1);frozen.XiX(i,j,k)=frozen.EtY(i,j,k)=1;
        }
        FDM::TurbulenceConfig tc;tc.enabled=true;tc.family=FDM::TurbulenceFamily::RAS;tc.model=FDM::TurbulenceModelKind::kOmegaSST;tc.laminarDynamicViscosity=.00001;
        Turbulence::Manager manager(tc,&boundary);require(manager.initialize(frozen),"SST wall service did not initialize.");
        auto& scalars=manager.scalarFields();
        for(auto& value:scalars.values(Turbulence::ScalarSlot::K))value=.01;
        for(auto& value:scalars.values(Turbulence::ScalarSlot::Omega))value=1;
        manager.correct(frozen,1e-8);
        const double cached=manager.dynamicViscosity(frozen,3,3,1,.00001);
        frozen(3,4,1,RU)=4;
        require(manager.dynamicViscosity(frozen,3,3,1,.00001)==cached,"Flow stage recomputed the frozen SST coefficient.");
        const int calls=boundary.calls;(void)manager.dynamicViscosity(frozen,3,3,1,.00001);
        require(boundary.calls==calls,"Viscosity view silently refreshed wall/transport state.");
        using namespace System;
        bool rejected=false;
        try { IBM::GhostIBM::applyViscousWall(frozen,{},1); }
        catch(const std::runtime_error&) {rejected=true;}
        require(rejected,"Invalid PerfectGas gamma was accepted by the wall closure.");
        rejected=false;
        try { IBM::GhostIBM::wallProjection({{1,0},{2,0}},{1,1}); }
        catch(const std::runtime_error&) {rejected=true;}
        require(rejected,"Rank-deficient wall support silently downgraded.");
        FDM::SolverConfig config;config.numerics.timeRecipe=FDM::builtInTimeRecipe(FDM::TimeRecipeId::ClassicalRK4);config.numerics.recipes.time=config.numerics.timeRecipe;
        BuildRequest request;request.singleFluidPreset=SingleFluidPresetSpec{true};request.composition.stateDeclared=true;request.composition.solutionVariables={"rho","rhoU","rhoE"};
        SystemContribution wall,sst;auto descriptor=IBM::Descriptor::ghostCell();descriptor.wallClosure=FDM::ImmersedWallClosure::StationaryNoSlipAdiabatic;
        IBM::SystemContribution::contribute(wall,descriptor);Turbulence::contribute(sst,{"kOmegaSST",{},false});request.modelContributions={wall,sst};
        const auto resolved=build(config,request);
        StateRegistry geometryState;geometryState.add(resolved.executableSystem.state.at("immersed.wallDistance"));
        State::StateBundle bundle;bundle.patches={&frozen};
        const auto geometryView=realizeState(geometryState,RuntimeRequirements{},bundle);
        const int cell=frozen.getIdx(3,3,1);
        require(geometryView.at("immersed.wallDistance").fields[0]->read(cell,0)==1,"Wall distance view does not alias geometry.");
        frozen.setWallDistance(3,3,1,2);
        require(geometryView.at("immersed.wallDistance").fields[0]->read(cell,0)==2,"Wall distance view copied storage.");
        for (const auto& binding:resolved.runtime.operationBindings)
            require(binding.status==BindingStatus::Resolved,"SST viscous immersed composition is not runnable.");
    }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
