#pragma once

/// @file SF_eulerian.h
/// @brief Eulerian 压力步进器内部结构网格算子。

#include "core/mesh/SF_dimension.h"
#include "SF_canonicalFace.h"
#include "SF_field.h"
#include "methods/numerics/structured/SF_structured.h"
#include "SF_scalarField.h"
#include "methods/numerics/structured/SF_vectorCalculus.h"

#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace SF::EulerianEulerian::Ops {

inline void offset(int axis,int sign,int& di,int& dj,int& dk){
    di=dj=dk=0;if(axis==0)di=sign;else if(axis==1)dj=sign;else dk=sign;
}

inline std::array<double,3> metric(const Field& f,int axis,int i,int j,int k){
    if(axis==0)return {f.XiX(i,j,k),f.XiY(i,j,k),f.XiZ(i,j,k)};
    if(axis==1)return {f.EtX(i,j,k),f.EtY(i,j,k),f.EtZ(i,j,k)};
    return {f.ZeX(i,j,k),f.ZeY(i,j,k),f.ZeZ(i,j,k)};
}

/// @brief 返回从 lower 点指向 upper 点的 canonical 面共因子面积向量。
inline std::array<double,3> faceCofactor(
        const Field& f,int axis,int lowerI,int lowerJ,int lowerK,
        bool axisymmetric=false,int radialCoordinate=0){
    std::array<double,3> result=
        Numerics::CanonicalFace::geometry(
            f,axis,lowerI,lowerJ,lowerK).cofactor;
    if(axisymmetric){
        int di=0,dj=0,dk=0;
        offset(axis,1,di,dj,dk);
        const auto coordinate=[&](int i,int j,int k){
            return radialCoordinate==0?f.X(i,j,k)
                :(radialCoordinate==1?f.Y(i,j,k):f.Z(i,j,k));
        };
        const double radius=0.5*(
            coordinate(lowerI,lowerJ,lowerK)
            +coordinate(lowerI+di,lowerJ+dj,lowerK+dk));
        if(!std::isfinite(radius)||radius<0.0)
            throw std::runtime_error(
                "Axisymmetric face has invalid radial coordinate: axis="
                +std::to_string(axis)+", lower=("
                +std::to_string(lowerI)+","+std::to_string(lowerJ)+","
                +std::to_string(lowerK)+"), radius="
                +std::to_string(radius)+".");
        constexpr double twoPi=6.283185307179586476925286766559;
        for(double& value:result)value*=twoPi*radius;
    }
    return result;
}

/// @brief 返回物理控制体积倒数 `1/V=Jac`，并拒绝退化网格。
inline double inverseCellVolume(
        const Field& f,int i,int j,int k,
        bool axisymmetric=false,int radialCoordinate=0){
    double inverseVolume=f.Jac(i,j,k);
    if(axisymmetric){
        const double radius=radialCoordinate==0?f.X(i,j,k)
            :(radialCoordinate==1?f.Y(i,j,k):f.Z(i,j,k));
        if(!std::isfinite(radius)||radius<=0.0)
            throw std::runtime_error(
                "Axisymmetric solved cell requires radius > 0.");
        constexpr double twoPi=6.283185307179586476925286766559;
        inverseVolume/=twoPi*radius;
    }
    if(!std::isfinite(inverseVolume)||inverseVolume<=0.0)
        throw std::runtime_error(
            "Eulerian operator found non-positive cell inverse volume.");
    return inverseVolume;
}

inline bool solved(
        const Field& f,int i,int j,int k,
        bool axisymmetric=false,int radialCoordinate=0){
    const int ng=f.NG();
    const bool base=i>=ng&&i<ng+f.NX()&&j>=ng&&j<ng+f.NY()
        &&k>=ng&&k<ng+f.NZ()&&f.CellFlag(i,j,k)==FLUID_CELL
        &&!f.isSolverBoundaryPoint(i,j,k);
    if(!base||!axisymmetric)return base;
    const double radius=radialCoordinate==0?f.X(i,j,k)
        :(radialCoordinate==1?f.Y(i,j,k):f.Z(i,j,k));
    if(!std::isfinite(radius)||radius<-1.0e-14)
        throw std::runtime_error(
            "Axisymmetric physical point has invalid negative radius.");
    return radius>1.0e-14;
}

inline double spacing(const Field& f,int i,int j,int k,
                      int ni,int nj,int nk){
    const double dx=f.X(ni,nj,nk)-f.X(i,j,k);
    const double dy=f.Y(ni,nj,nk)-f.Y(i,j,k);
    const double dz=f.Z(ni,nj,nk)-f.Z(i,j,k);
    const double h=std::sqrt(dx*dx+dy*dy+dz*dz);
    if(!std::isfinite(h)||h<=0.0)throw std::runtime_error("Eulerian operator found degenerate spacing.");
    return h;
}

inline Vector3 gradient(const Field& f,const ScalarField& q,
                        int i,int j,int k){
    return CENTRAL2::grad(f, q, i, j, k);
}

inline double divergenceFlux(const Field& f,const ScalarField& mass,
                             const std::array<ScalarField,3>& velocity,
                             int i,int j,int k){
    double div=0.0;
    for(int a=0;a<3;++a){
        if(!Math::isDirectionActiveIndex(a))continue;
        int ip,jp,kp,im,jm,km;offset(a,1,ip,jp,kp);offset(a,-1,im,jm,km);
        const auto g=metric(f,a,i,j,k);
        for(int d=0;d<3;++d){
            const double fp=mass(i+ip,j+jp,k+kp)*velocity[(size_t)d](i+ip,j+jp,k+kp);
            const double fm=mass(i+im,j+jm,k+km)*velocity[(size_t)d](i+im,j+jm,k+km);
            div+=0.5*g[(size_t)d]*(fp-fm);
        }
    }
    return div;
}

} // namespace SF::EulerianEulerian::Ops
