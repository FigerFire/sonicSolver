#pragma once

/// @file SF_nodalQuadrature.h
/// @brief source-canonical 结构网格点的物理对偶体积，不把分区端点当成物理端点。
#include "core/field/SF_field.h"
#include <cmath>
#include <stdexcept>

namespace SF::StructuredMesh {
/// @brief 梯形积分的 nodal volume；EMPTY 两端仍共同表示同一物理厚度。
inline double nodalVolume(const Field& field,int i,int j,int k) {
    const double jac=field.Jac(i,j,k);
    if(!std::isfinite(jac)||jac==0.) throw std::runtime_error("Nodal quadrature requires a finite nonzero Jacobian.");
    double weight=1.;
    const int lo=field.NG();
    const int coordinate[3]={i,j,k};const int size[3]={field.NX(),field.NY(),field.NZ()};
    for(int axis=0;axis<3;++axis) {
        if(size[axis]<=1)continue;
        int offset=0;
        if(coordinate[axis]==lo)offset=-1;
        else if(coordinate[axis]==lo+size[axis]-1)offset=1;
        if(!offset)continue;
        const int ni=i+(axis==0?offset:0),nj=j+(axis==1?offset:0),nk=k+(axis==2?offset:0);
        // Interface halo metadata comes from the source mesh topology. Shared
        // Eulerian replicas never acquire an extra half-volume due to partitioning.
        if(lo==0 || !field.isCommunicationHalo(ni,nj,nk))weight*=.5;
    }
    return weight/std::abs(jac);
}
}
