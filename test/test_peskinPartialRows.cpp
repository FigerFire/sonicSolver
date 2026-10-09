#include "models/ibm/method/SF_method.h"
#include "infrastructure/mpi/SF_parallelContext.h"
#include "infrastructure/execution/SF_executionRuntime.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
class PartialSurface final : public SF::IBM::Forcing::ISurfaceConstraintOperator {
public:
    PartialSurface(bool contributes,int cell):contributes_(contributes),cell_(cell){}
    const SF::FDM::ImmersedSurfaceSystem& build(const SF::Field&,
        const std::vector<SF::IBM::GeoProcessing::Triangle>&,const SF::Vector3&,double,SF::FDM::IBMSurfaceNormalization) override {
        normalization_={contributes_?1.:0.};return system_;
    }
    const std::vector<double>& localNormalizations() const override {return normalization_;}
    const SF::FDM::ImmersedSurfaceSystem& normalizeDistributed(const std::vector<double>& normalization,const std::vector<double>&) override {
        if(normalization.size()!=1 || normalization[0]!=1) throw std::runtime_error("partial graph normalization failed");
        system_.points.clear();SF::FDM::ImmersedSurfacePoint marker;
        marker.globalConstraintId=SF::GlobalConstraintDofId::fromMarker(SF::GlobalMarkerId::fromSurfacePrimitive(7));
        marker.measure=.7;
        if(contributes_) {SF::FDM::ImmersedInterpolationWeight edge;edge.cell=cell_;edge.value=1.;edge.dualVolume=3.;edge.globalEulerianDofId=0;marker.interpolation.push_back(edge);}
        system_.points.push_back(marker);return system_;
    }
private:
    bool contributes_;int cell_;std::vector<double> normalization_;SF::FDM::ImmersedSurfaceSystem system_;
};
}
int main(int argc,char** argv) {
    try {
        SF::Parallel::ParallelContext parallel(argc,argv,true);
        SF::Execution::Runtime runtime(&parallel.coordinator());
        SF::Field field;field.setup(3,3,1,1,5);const int cell=field.getIdx(2,2,1);
        field(2,2,1,SF::RHO)=2.;field(2,2,1,SF::RU)=.8;field(2,2,1,SF::E)=10.;field.Jac(2,2,1)=1./3.;
        SF::IBM::GeoProcessing::STLGeometry geometry;
        if(argc<2 || !geometry.loadSTLFiles({argv[1]},".")) throw std::runtime_error("test surface STL unavailable");
        SF::IBM::IBMRuntimeConfig config;auto& forcing=config.forcing;
        forcing.algorithm=SF::FDM::IBMForcingAlgorithm::PeskinOriginal;
        forcing.constraintSupport=SF::FDM::IBMConstraintSupport::Surface;
        forcing.constraintDomain=SF::FDM::IBMConstraintDomain::Surface;
        forcing.representation=SF::FDM::IBMRepresentation::DiffuseKernel;
        forcing.enforcement=SF::FDM::IBMEnforcement::ExplicitIBM;
        forcing.surfaceKernel="wendlandC2";forcing.surfaceQuadrature="triangleCentroid";
        forcing.surfaceNormalization=SF::FDM::IBMSurfaceNormalization::PartitionOfUnity;forcing.surfaceSpreading="adjoint";forcing.surfaceSupportRadius=2.;
        SF::IBM::Forcing::ImmersedForcingSystem method;
        method.configure(geometry,config);method.setExecutionRuntime(&runtime);
        method.setSurfaceConstraintOperator(std::make_unique<PartialSurface>(parallel.rank()==0,cell));
        const auto first=method.projectPredictedState({&field},.1,.1);
        if(first.constrainedCells!=1 || field(2,2,1,SF::RU)!=.8 || first.maximumVelocityResidual!=.4)
            throw std::runtime_error("Peskin ceased to be a lagged canonical marker response");
        const auto second=method.projectPredictedState({&field},.2,.1);
        if(parallel.rank()==0 && (std::abs(field(2,2,1,SF::RU))>1e-14 || std::abs(field(2,2,1,SF::E)-9.84)>1e-13))
            throw std::runtime_error("lagged force or midpoint energy not applied to Eulerian owner");
        if(std::abs(second.forceOnBody.x-24.)>1e-12 || std::abs(second.fluidMechanicalPower+4.8)>1e-12)
            throw std::runtime_error("Peskin owner load or volume-weighted power duplicated");
        if(parallel.isRoot())std::cout<<"actual Peskin empty partial rows, remote marker owner, lagged update and midpoint power passed\n";
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
