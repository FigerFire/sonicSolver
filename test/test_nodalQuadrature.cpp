#include "core/mesh/SF_nodalQuadrature.h"
#include <iostream>
#include <stdexcept>
int main() {
    try {
        SF::Field field;field.setup(5,3,2,1,5);
        for(int c=0;c<field.TotalSize();++c) {int i,j,k;field.getIJK(c,i,j,k);field.Jac(i,j,k)=10.;}
        using SF::StructuredMesh::nodalVolume;
        if(nodalVolume(field,3,2,1)!=.05 || nodalVolume(field,5,2,1)!=.025 || nodalVolume(field,1,1,1)!=.0125)
            throw std::runtime_error("physical endpoint/EMPTY volume is incorrect");
        field.setCommunicationHalo(6,2,1,true);
        if(nodalVolume(field,5,2,1)!=.05)throw std::runtime_error("partition endpoint incorrectly halved volume");
        double v=0.;
        field.clearCommunicationHaloMask();
        for(int k=1;k<=2;++k)for(int j=1;j<=3;++j)for(int i=1;i<=5;++i)v+=nodalVolume(field,i,j,k);
        if(std::abs(v-.8)>1e-14)throw std::runtime_error("integrated nodal volume differs from physical domain volume");
        SF::Field noHalo;noHalo.setup(3,3,2,0,5);
        noHalo.Jac(0,0,0)=10.;
        if(nodalVolume(noHalo,0,0,0)!=.0125)
            throw std::runtime_error("physical corner quadrature requires no halo storage");
        std::cout<<"physical nodal and partition endpoint quadrature passed\n";
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
