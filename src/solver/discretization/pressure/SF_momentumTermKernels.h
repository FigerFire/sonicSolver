#pragma once

/// @file SF_momentumTermKernels.h
/// @brief Existing primitive Upwind1/Central2 pressure-Momentum stencils.

#include "core/field/SF_field.h"
#include "core/state/SF_distributedField.h"
#include "core/mesh/SF_dimension.h"
#include "methods/numerics/structured/SF_canonicalFace.h"

#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace SF::Discretization::PressureMomentum {

struct Stencil {
    const State::DistributedFieldView& velocity;
    const Field& geometry;
    const std::vector<double>& faceFlux;
};
using Kernel = double (*)(const Stencil&, int cell, int component);

namespace detail {

inline int adjacent(const Field& geometry,int cell,int axis,int sign) {
    int i=0,j=0,k=0;
    geometry.getIJK(cell,i,j,k);
    if (axis==0) i+=sign;
    else if (axis==1) j+=sign;
    else k+=sign;
    if (i<0 || i>=geometry.MX() || j<0 || j>=geometry.MY()
        || k<0 || k>=geometry.MZ())
        throw std::runtime_error("Pressure stencil exceeds allocated halo.");
    return geometry.getIdx(i,j,k);
}

inline std::size_t faceIndex(const Field& geometry,int axis,int lower) {
    return (size_t)axis*(size_t)geometry.TotalSize()+(size_t)lower;
}

inline std::array<double,3> faceArea(const Field& geometry,
                                      int axis,int lower) {
    int i=0,j=0,k=0;
    geometry.getIJK(lower,i,j,k);
    return Numerics::CanonicalFace::geometry(geometry,axis,i,j,k).cofactor;
}

inline double spacing(const Field& geometry,int axis,int lower) {
    const int upper=adjacent(geometry,lower,axis,1);
    int i=0,j=0,k=0,ni=0,nj=0,nk=0;
    geometry.getIJK(lower,i,j,k);
    geometry.getIJK(upper,ni,nj,nk);
    const double dx=geometry.X(ni,nj,nk)-geometry.X(i,j,k);
    const double dy=geometry.Y(ni,nj,nk)-geometry.Y(i,j,k);
    const double dz=geometry.Z(ni,nj,nk)-geometry.Z(i,j,k);
    const double h=std::sqrt(dx*dx+dy*dy+dz*dz);
    if (!std::isfinite(h) || h<=0.0)
        throw std::runtime_error("Pressure operator found degenerate face spacing.");
    return h;
}

} // namespace detail

inline double primitiveUpwind1(const Stencil& stencil,
                               int cell,int component) {
    double advective=0.0;
    for (int axis=0;axis<3;++axis) {
        if (!Math::isDirectionActiveIndex(axis)) continue;
        for (int sign:{-1,1}) {
            const int lower=sign>0?cell:detail::adjacent(stencil.geometry,cell,axis,-1);
            const int upper=detail::adjacent(stencil.geometry,lower,axis,1);
            const double flux=stencil.faceFlux[
                detail::faceIndex(stencil.geometry,axis,lower)];
            const double upwind=stencil.velocity.read(
                flux>=0.0?lower:upper,component);
            advective+=sign*flux*upwind;
        }
    }
    return advective;
}

inline double central2(const Stencil& stencil,
                       int cell,int component) {
    double laplacian=0.0;
    for (int axis=0;axis<3;++axis) {
        if (!Math::isDirectionActiveIndex(axis)) continue;
        for (int sign:{-1,1}) {
            const int lower=sign>0?cell:detail::adjacent(stencil.geometry,cell,axis,-1);
            const auto area=detail::faceArea(stencil.geometry,axis,lower);
            const double magnitude=std::sqrt(
                area[0]*area[0]+area[1]*area[1]+area[2]*area[2]);
            const int neighbour=detail::adjacent(stencil.geometry,cell,axis,sign);
            laplacian+=magnitude
                *(stencil.velocity.read(neighbour,component)
                  -stencil.velocity.read(cell,component))
                /detail::spacing(stencil.geometry,axis,lower);
        }
    }
    return laplacian;
}

} // namespace SF::Discretization::PressureMomentum
