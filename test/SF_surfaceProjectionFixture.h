#pragma once
#include "models/ibm/constraint/variational/SF_projection.h"
#include "methods/numerics/immersed/SF_immersed.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <functional>

inline void verifySurfaceProjection(int rank=0,int ranks=1,
    std::function<void(std::vector<double>&)> sum={}) {
    using namespace SF;
    FDM::ImmersedSurfaceSystem surface;
    const double weights[2][3]={{.75,.25,0.},{0.,.5,.5}};
    const double volume[3]={.3,.7,1.1},rho[3]={1.2,.9,1.4};
    const Vector3 u[3]={{.2,.4,-.1},{.8,-.1,.3},{-.3,.7,.5}};
    const Vector3 lambda[2]={{1.,-.5,.2},{2.,.4,-.6}};
    std::vector<double> inverseMass(3,0.0);
    for(int cell=0;cell<3;++cell) if(cell%ranks==rank) inverseMass[cell]=1/(rho[cell]*volume[cell]);
    for(int marker=0;marker<2;++marker) {
        FDM::ImmersedSurfacePoint point;point.measure=marker==0?.4:.9;
        point.globalConstraintId=GlobalConstraintDofId::fromMarker(GlobalMarkerId::fromSurfacePrimitive(marker));
        for(int cell=0;cell<3;++cell) if(weights[marker][cell]>0 && cell%ranks==rank) {
            FDM::ImmersedInterpolationWeight edge;edge.cell=cell;edge.value=weights[marker][cell];edge.globalEulerianDofId=cell;edge.dualVolume=volume[cell];point.interpolation.push_back(edge);
        }
        surface.points.push_back(point);
    }
    std::vector<double> unity(2,0.0),ju(6,0.0),diagonal(2,0.0);
    std::vector<Vector3> force(3);
    for(int marker=0;marker<2;++marker) {
        const auto& point=surface.points[marker];
        const Vector3 value=FDM::Immersed::interpolate(point.interpolation,[&](int cell){return u[cell];});
        ju[3*marker]=value.x;ju[3*marker+1]=value.y;ju[3*marker+2]=value.z;
        for(auto edge:point.interpolation) {unity[marker]+=edge.value;diagonal[marker]+=edge.value*edge.value*point.measure*inverseMass[edge.cell];}
        FDM::Immersed::spread(point.interpolation,point.measure,lambda[marker],[&](int cell,Vector3 f){force[cell]=force[cell]+f*(1/volume[cell]);});
    }
    // Ranks outside support legally contribute all-zero partial diagonals.
    for(auto value:diagonal) if(!std::isfinite(value)||value<0) throw std::runtime_error("invalid partial diagonal");
    if(sum) {sum(unity);sum(ju);sum(diagonal);}
    double adjointL=0,adjointR=0;
    for(int marker=0;marker<2;++marker) {
        if(std::abs(unity[marker]-1)>1e-14 || diagonal[marker]<=0) throw std::runtime_error("global unity or positive mass response failed");
        if(marker%ranks==rank) adjointL+=surface.points[marker].measure*dot(Vector3{ju[3*marker],ju[3*marker+1],ju[3*marker+2]},lambda[marker]);
    }
    for(int cell=0;cell<3;++cell) if(cell%ranks==rank) adjointR+=volume[cell]*dot(u[cell],force[cell]);
    std::vector<double> adjoint{adjointL,adjointR};if(sum)sum(adjoint);
    if(std::abs(adjoint[0]-adjoint[1])>1e-14) throw std::runtime_error("weighted J/Jtranspose adjoint failed");
    // Form rhs from a known exact multiplier using the complete coupled graph.
    IBM::Variational::SurfaceSchurProjectionProblem problem;problem.surface=&surface;problem.inverseEulerianMass=&inverseMass;
    problem.dt=.1;problem.maxIterations=20;problem.relativeTolerance=1e-12;problem.globalSum=sum;
    for(int marker=0;marker<2;++marker) problem.ownedMarkers.push_back(marker%ranks==rank);
    std::vector<double> rhs(6,0.0);
    for(int marker=0;marker<2;++marker) for(auto edge:surface.points[marker].interpolation) {
        auto corrected=force[edge.cell]*(problem.dt/rho[edge.cell]);
        rhs[3*marker]+=edge.value*corrected.x;rhs[3*marker+1]+=edge.value*corrected.y;rhs[3*marker+2]+=edge.value*corrected.z;
    }
    if(sum)sum(rhs);
    for(int marker=0;marker<2;++marker) problem.constraintError.push_back({rhs[3*marker],rhs[3*marker+1],rhs[3*marker+2]});
    const auto solution=IBM::Variational::solveSurfaceSchurProjection(problem);
    for(int marker=0;marker<2;++marker) if(norm(solution.multiplier[marker]-lambda[marker])>1e-10) throw std::runtime_error("coupled Schur solution differs from exact lambda");
    if(rank==0) std::cout<<"surface unity, weighted adjoint, partial diagonal and Schur projection passed\n";
}
