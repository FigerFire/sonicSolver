#include "solver/boundary/reconstruction/SF_reconstruction.h"
#include "solver/boundary/reconstruction/ILW/SF_boundaryClosure.h"
#include "core/mesh/SF_dimension.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void initialize(SF::Field& field,int nx,int ny) {
    field.setup(nx,ny,2,3,5);
    for(int k=0;k<field.MZ();++k)for(int j=0;j<field.MY();++j)for(int i=0;i<field.MX();++i) {
        const double x=i-field.NG(),y=j-field.NG();
        field.X(i,j,k)=x;field.Y(i,j,k)=y;field.Z(i,j,k)=.1*(k-field.NG());
        field.Jac(i,j,k)=10.;
        field(i,j,k,0)=1.2+.01*std::sin(.3*x)*std::cos(.2*y);
    }
}
}
int main() {
    try {
        SF::Math::resetActiveDirections();SF::Math::deactivateDirection(2);
        for(int order:{3,5}) for(int cutAxis:{0,1}) {
            SF::Field serial,partition;
            initialize(serial,cutAxis==0?21:9,cutAxis==0?9:21);
            initialize(partition,cutAxis==0?11:9,cutAxis==0?9:11);
            const int ng=partition.NG(),interface=ng+10,wallAxis=1-cutAxis;
            for(int k=ng;k<ng+partition.NZ();++k)
                for(int tangent=ng;tangent<ng+9;++tangent)
                    for(int layer=1;layer<=ng;++layer)
                        partition.setCommunicationHalo(cutAxis==0?interface+layer:tangent,
                            cutAxis==1?interface+layer:tangent,k,true);
            const int i=cutAxis==0?interface:ng,j=cutAxis==1?interface:ng;
            const auto normal=SF::Boundary::ILW::BoundaryClosure::boundaryNormal(partition,i,j,ng);
            if(normal.axis!=wallAxis || normal.sign!=-1)
                throw std::runtime_error("partition endpoint was treated as the physical wall normal");
            for(int layer=1;layer<=ng;++layer) {
                const int gi=i-(wallAxis==0?layer:0),gj=j-(wallAxis==1?layer:0);
                serial(gi,gj,ng,0)=-999.;partition(gi,gj,ng,0)=-999.;
            }
            SF::Boundary::Reconstruction::Selection recipe{SF::Boundary::Reconstruction::Kind::ILW,order};
            for(auto* field:{&serial,&partition})
                SF::Boundary::Reconstruction::applyScalar(*field,i,j,ng,wallAxis,0,0.,
                    SF::Boundary::Law::Kind::ZeroGradient,recipe);
            for(int layer=1;layer<=ng;++layer) {
                const int gi=i-(wallAxis==0?layer:0),gj=j-(wallAxis==1?layer:0);
                const double value=serial(gi,gj,ng,0);
                if(!std::isfinite(value) || value<=0. || value!=partition(gi,gj,ng,0))
                    throw std::runtime_error("physical ILW ghost differs from the same unpartitioned stencil");
            }
        }
        std::cout<<"physical normal and order 3/5 tangential ILW stencils are partition invariant\n";
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
