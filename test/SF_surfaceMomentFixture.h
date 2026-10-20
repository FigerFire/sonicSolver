#pragma once
#include "models/ibm/topology/SF_couplingGraph.h"
#include "models/ibm/constraint/variational/SF_projection.h"
#include "methods/numerics/immersed/SF_immersed.h"
#include "core/mesh/SF_nodalQuadrature.h"
#include <algorithm>
#include <functional>
#include <iostream>

inline void verifySurfaceMoments(int rank=0,int ranks=1,
        std::function<void(std::vector<double>&)> sum={}) {
    using namespace SF;
    using N=FDM::IBMSurfaceNormalization;
    auto require=[](bool valid,const char* message) {if (!valid) throw std::runtime_error(message);};
    Field field;field.setup(7,7,2,1,5);field.clearCellFlags(FLUID_CELL);
    for (int k=1;k<=2;++k) for (int j=1;j<=7;++j) for (int i=1;i<=7;++i) {
        field.X(i,j,k)=.4*(i-1);field.Y(i,j,k)=.3*(j-1);field.Z(i,j,k)=.1*(k-1);
        field.Jac(i,j,k)=1./(.012*(1.+.02*i+.01*j));
        const int cell=field.getIdx(i,j,k),owner=cell%std::min(ranks,2);
        field.setGlobalDofOwnership(i,j,k,cell,owner,owner==rank);
        // One-sided wall support excludes the physical x=0 node, as the real
        // graph excludes boundary unknowns. Markers outside its convex hull
        // require signed weights; clipping would invalidate first moments.
        if (i==1) field.CellFlag(i,j,k)=SOLID_CELL;
    }
    bool sawNegative=false;
    double maxMoment=0.,maxTorque=0.,maxPower=0.,maxPartition=0.;
    for (double shift:{0.,.013,.17,.38,.8}) {
        std::vector<IBM::Topology::SurfaceMarker> markers(2);
        markers[0].id=GlobalMarkerId::fromSurfacePrimitive(0);markers[0].position={.1+shift,.91,.043};markers[0].measure=.3;
        markers[1].id=GlobalMarkerId::fromSurfacePrimitive(1);markers[1].position={.73+shift,1.15,.06};markers[1].measure=.9;
        IBM::Topology::CouplingGraph graph;graph.build(field,markers,1.3,N::LinearReproducing);
        auto normalization=graph.rawNormalizations(),moments=graph.rawMoments();
        if(sum) {sum(normalization);sum(moments);}graph.normalize(normalization,moments);
        const auto& rows=graph.rows();const int total=field.TotalSize();
        auto position=[&](int cell) {int i,j,k;field.getIJK(cell,i,j,k);return Vector3{field.X(i,j,k),field.Y(i,j,k),field.Z(i,j,k)};};
        const Vector3 center{.65,.52,.03},omega{.3,-.5,1.4};
        auto rotation=[&](Vector3 x) {return Vector3{.2,-.3,.1}+cross(omega,x-center);};
        auto affine=[&](Vector3 x) {return Vector3{1.+2*x.x-x.y+.3*x.z,.7-.2*x.x+.4*x.y-x.z,-.2+x.x+.2*x.z};};
        auto velocity=[&](int cell) {const auto x=position(cell);return Vector3{std::sin(x.x)+x.y*x.y,x.z+x.x*x.y,std::cos(x.y)+x.z*x.z};};
        const Vector3 lambda[2]={{.7,-1.3,.5},{-.2,.4,.6}};
        std::vector<Vector3> force(total);std::vector<double> inverseMass(total,0.);
        std::vector<double> reproduced(2*10,0.),weights(2*total,0.);
        FDM::ImmersedSurfaceSystem surface;
        for (int m=0;m<2;++m) {
            FDM::ImmersedSurfacePoint point;point.position=markers[m].position;point.measure=markers[m].measure;point.interpolation=rows[m];surface.points.push_back(point);
            for (const auto& edge:rows[m]) {
                sawNegative=sawNegative || edge.value<0;
                const auto x=position(edge.cell),a=affine(x),r=rotation(x);
                reproduced[10*m]+=edge.value;
                for(int c=0;c<3;++c) {
                    const double px=c==0?x.x:c==1?x.y:x.z;
                    const double av=c==0?a.x:c==1?a.y:a.z,rv=c==0?r.x:c==1?r.y:r.z;
                    reproduced[10*m+1+c]+=edge.value*px;
                    reproduced[10*m+4+c]+=edge.value*av;
                    reproduced[10*m+7+c]+=edge.value*rv;
                }
                inverseMass[edge.cell]=1./((.9+.001*edge.cell)*edge.dualVolume);
                weights[m*total+edge.cell]=edge.value;
            }
            FDM::Immersed::spread(rows[m],markers[m].measure,lambda[m],
                [&](int cell,Vector3 value) {int i,j,k;field.getIJK(cell,i,j,k);force[cell]=force[cell]+value*(1./StructuredMesh::nodalVolume(field,i,j,k));});
        }
        if(sum) {sum(reproduced);sum(weights);}
        for (int m=0;m<2;++m) {
            const auto x=markers[m].position,a=affine(x),r=rotation(x);
            require(std::abs(reproduced[10*m]-1.)<3e-12,"constant reproduction failed");
            for(int c=0;c<3;++c) {
                const double px=c==0?x.x:c==1?x.y:x.z,av=c==0?a.x:c==1?a.y:a.z,rv=c==0?r.x:c==1?r.y:r.z;
                maxMoment=std::max(maxMoment,std::abs(reproduced[10*m+1+c]-px));
                require(std::abs(reproduced[10*m+4+c]-av)<3e-12,"linear field reproduction failed");
                require(std::abs(reproduced[10*m+7+c]-rv)<3e-12,"rigid rotation interpolation failed");
            }
        }
        // Independent serial assembly from all unique nodes checks that global
        // moments, not partition-local fits, select the distributed weights.
        Field full;full.setup(7,7,2,1,5);full.clearCellFlags(FLUID_CELL);
        for (int k=1;k<=2;++k) for(int j=1;j<=7;++j) for(int i=1;i<=7;++i) {
            full.X(i,j,k)=field.X(i,j,k);full.Y(i,j,k)=field.Y(i,j,k);full.Z(i,j,k)=field.Z(i,j,k);full.Jac(i,j,k)=field.Jac(i,j,k);
            if(i==1)full.CellFlag(i,j,k)=SOLID_CELL;
        }
        IBM::Topology::CouplingGraph serial;serial.build(full,markers,1.3,N::LinearReproducing);
        serial.normalize(serial.rawNormalizations(),serial.rawMoments());
        for(int m=0;m<2;++m) for(const auto& edge:serial.rows()[m]) maxPartition=std::max(maxPartition,std::abs(edge.value-weights[m*total+edge.cell]));
        std::vector<double> balance(16,0.); // Eulerian force/torque/power and marker counterparts.
        for(int cell=0;cell<total;++cell) if(inverseMass[cell]>0.) {
            int i,j,k;field.getIJK(cell,i,j,k);const double volume=StructuredMesh::nodalVolume(field,i,j,k);
            const auto impulse=force[cell]*volume,moment=cross(position(cell)-center,impulse),u=velocity(cell);
            balance[0]+=impulse.x;balance[1]+=impulse.y;balance[2]+=impulse.z;
            balance[3]+=moment.x;balance[4]+=moment.y;balance[5]+=moment.z;balance[6]+=dot(u,impulse);
            // For arbitrary dt/rho, the same impulse is the angular momentum
            // exchange and midpoint work equals the kinetic-energy increment.
            const double dt=.017,rho=.9+.001*cell;
            const auto updated=u+force[cell]*(dt/rho);
            balance[14]+=.5*rho*volume*(dot(updated,updated)-dot(u,u));
            balance[15]+=dt*dot((u+updated)*.5,impulse);
        }
        std::vector<double> ju(6,0.);
        for(int m=0;m<2;++m) {
            const auto value=FDM::Immersed::interpolate(rows[m],velocity);ju[3*m]=value.x;ju[3*m+1]=value.y;ju[3*m+2]=value.z;
        }
        if(sum)sum(ju);
        for(int m=0;m<2;++m) if(m%ranks==rank) {
            const auto load=lambda[m]*markers[m].measure,moment=cross(markers[m].position-center,load);
            balance[7]+=load.x;balance[8]+=load.y;balance[9]+=load.z;
            balance[10]+=moment.x;balance[11]+=moment.y;balance[12]+=moment.z;
            balance[13]+=dot(Vector3{ju[3*m],ju[3*m+1],ju[3*m+2]},load);
        }
        if(sum)sum(balance);
        for(int c=0;c<3;++c) {
            require(std::abs(balance[c]-balance[7+c])<3e-12,"force exchange failed");
            maxTorque=std::max(maxTorque,std::abs(balance[3+c]-balance[10+c]));
        }
        maxPower=std::max(maxPower,std::abs(balance[6]-balance[13]));
        require(std::abs(balance[14]-balance[15])<3e-12,"mechanical work / kinetic increment failed");
        IBM::Variational::SurfaceSchurProjectionProblem problem;
        problem.surface=&surface;problem.inverseEulerianMass=&inverseMass;problem.dt=.017;
        problem.maxIterations=30;problem.relativeTolerance=1e-12;problem.globalSum=sum;
        std::vector<double> rhs(6,0.);
        for(int m=0;m<2;++m) {
            const auto value=FDM::Immersed::interpolate(rows[m],[&](int cell) {return force[cell]*(problem.dt/(.9+.001*cell));});
            rhs[3*m]=value.x;rhs[3*m+1]=value.y;rhs[3*m+2]=value.z;problem.ownedMarkers.push_back(m%ranks==rank);
        }
        if(sum)sum(rhs);
        for(int m=0;m<2;++m)problem.constraintError.push_back({rhs[3*m],rhs[3*m+1],rhs[3*m+2]});
        const auto solved=IBM::Variational::solveSurfaceSchurProjection(problem);
        for(int m=0;m<2;++m)require(norm(solved.multiplier[m]-lambda[m])<1e-10,"signed-weight Schur solve failed");
    }
    require(maxMoment<3e-12 && maxTorque<3e-12 && maxPower<3e-12 && maxPartition<3e-12,"transfer moment/torque/adjoint/partition contract failed");
    std::vector<double> negative{double(sawNegative)};if(sum)sum(negative);
    require(negative[0]>0,"one-sided test never exercised signed weights");
    // A planar support cannot represent off-plane affine variations. Reject it
    // on every rank after complete moment SUM instead of using a lower-order fit.
    for(int k=1;k<=2;++k)for(int j=1;j<=7;++j)for(int i=1;i<=7;++i)field.Z(i,j,k)=.04;
    std::vector<IBM::Topology::SurfaceMarker> marker(1);marker[0].position={1.,.8,.04};
    IBM::Topology::CouplingGraph planar;planar.build(field,marker,1.3,N::LinearReproducing);
    auto norm=planar.rawNormalizations(),mom=planar.rawMoments();if(sum){sum(norm);sum(mom);}
    bool rejected=false;try{planar.normalize(norm,mom);}catch(const std::runtime_error&){rejected=true;}
    require(rejected,"rank-deficient support silently fell back");
    if(rank==0) std::cout<<"affine moment="<<maxMoment<<" torque="<<maxTorque<<" weighted adjoint="<<maxPower<<" partition="<<maxPartition<<"; signed weights, translation and rank rejection passed\n";
}
