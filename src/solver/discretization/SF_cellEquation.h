#pragma once

/// @file SF_cellEquation.h
/// @brief Cell-vector Upwind1/Central2 contributions and explicit numerical update.

#include "core/field/SF_field.h"
#include "core/state/SF_distributedField.h"
#include "core/mesh/SF_dimension.h"
#include "core/interfaces/SF_termKernel.h"
#include "methods/numerics/structured/SF_vectorCalculus.h"
#include "methods/math/SF_newtonian.h"
#include "methods/numerics/structured/SF_canonicalFace.h"

#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace SF::Discretization::CellEquation {

struct Stencil {
    const State::DistributedFieldView& velocity;
    const Field& geometry;
    const std::vector<double>& faceFlux;
    std::array<bool,3> active{true,true,true};
    const State::DistributedFieldView* coefficient=nullptr;
    double constantCoefficient=1;
    const std::vector<double>* traction=nullptr;
};
using Kernel = double (*)(const Stencil&, int cell, int component);
using GradientKernel=std::array<double,3> (*)(const Field&,const State::DistributedFieldView&,int);
inline std::array<double,3> centralGradient(const Field& geometry,const State::DistributedFieldView& value,int cell) {
    int i,j,k;geometry.getIJK(cell,i,j,k);
    const auto g=CENTRAL2::gradSampled(geometry,[&](int x,int y,int z){return value.read(geometry.getIdx(x,y,z),0);},i,j,k);
    return {g.x,g.y,g.z};
}
struct Contribution {Kernel spatial=nullptr;System::PrimitiveMomentumSource point;double weight=1;const State::DistributedFieldView* coefficient=nullptr;double constantCoefficient=1;};
struct GradientContribution {
    const State::DistributedFieldView* input;
    GradientKernel kernel;double weight,divisor;
};

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
        if (!stencil.active[axis]) continue;
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
        if (!stencil.active[axis]) continue;
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

/// @brief 面牵引的守恒散度；算术面系数和切向梯度平均是明确所选数值方法。
inline Vector3 stressFace(const Stencil& stencil,int axis,int lower) {
    const int upper=detail::adjacent(stencil.geometry,lower,axis,1);
    const auto mu=[&](int c){return stencil.coefficient?stencil.coefficient->read(c,0):stencil.constantCoefficient;};
    const double left=mu(lower),right=mu(upper);
    if(!std::isfinite(left) || left<0 || !std::isfinite(right) || right<0)throw std::runtime_error("Stress coefficient must be finite and nonnegative");
    const double faceMu=.5*(left+right);
    int i,j,k,ni,nj,nk;stencil.geometry.getIJK(lower,i,j,k);stencil.geometry.getIJK(upper,ni,nj,nk);
    Tensor3 gradient;
    const auto area=detail::faceArea(stencil.geometry,axis,lower);
    const double magnitude=std::sqrt(area[0]*area[0]+area[1]*area[1]+area[2]*area[2]);
    if(!std::isfinite(magnitude) || magnitude<=0)throw std::runtime_error("Stress face needs positive finite geometry area");
    const double h=detail::spacing(stencil.geometry,axis,lower);
    const Vector3 normal((stencil.geometry.X(ni,nj,nk)-stencil.geometry.X(i,j,k))/h,
        (stencil.geometry.Y(ni,nj,nk)-stencil.geometry.Y(i,j,k))/h,
        (stencil.geometry.Z(ni,nj,nk)-stencil.geometry.Z(i,j,k))/h);
    for(int c=0;c<3;++c) {
        const auto sample=[&](int x,int y,int z){return stencil.velocity.read(stencil.geometry.getIdx(x,y,z),c);};
        gradient[c]=.5*(CENTRAL2::gradSampled(stencil.geometry,sample,i,j,k)+CENTRAL2::gradSampled(stencil.geometry,sample,ni,nj,nk));
        // Replace only the interpolated normal derivative by the two-point face derivative.
        gradient[c]=gradient[c]+normal*((stencil.velocity.read(upper,c)-stencil.velocity.read(lower,c))/h-dot(gradient[c],normal));
    }
    return Constitutive::Newtonian::traction(Constitutive::Newtonian::stress(faceMu,gradient),Vector3(area[0],area[1],area[2]));
}
inline double stressArithmeticCentral2(const Stencil& stencil,int cell,int component) {
    if(!stencil.traction)throw std::runtime_error("Stress kernel lacks its once-evaluated face traction workspace");
    double divergence=0;
    for(int axis=0;axis<3;++axis)if(stencil.active[axis])for(int sign:{-1,1}) {
        const int lower=sign>0?cell:detail::adjacent(stencil.geometry,cell,axis,-1);
        divergence+=sign*stencil.traction->at((size_t)(axis*3+component)*stencil.geometry.TotalSize()+lower);
    }
    return divergence;
}

/// @brief 完整空间 RHS；只读所有邻域，不乘 dt、不写未知量。
inline void evaluate(const Stencil& stencil,const std::vector<int>& cells,const std::vector<Contribution>& terms,
        const std::vector<GradientContribution>& gradients,std::vector<double>& rhs,std::vector<double>* withoutGradient=nullptr) {
    const auto& field=stencil.geometry;const int total=field.TotalSize();rhs.assign((size_t)total*3,0);
    if(withoutGradient)withoutGradient->assign(rhs.size(),0);
    std::vector<std::vector<double>> traction(terms.size());
    for(std::size_t n=0;n<terms.size();++n)if(terms[n].spatial==&stressArithmeticCentral2) {
        auto bound=stencil;bound.coefficient=terms[n].coefficient;bound.constantCoefficient=terms[n].constantCoefficient;
        traction[n].resize((size_t)total*9);std::vector<bool> evaluated((size_t)total*3,false);
        // Each oriented face has one traction; adjacent cells consume opposite signs.
        for(int cell:cells)for(int axis=0;axis<3;++axis)if(stencil.active[axis])for(int sign:{-1,1}) {
            const int lower=sign>0?cell:detail::adjacent(field,cell,axis,-1);
            const auto face=detail::faceIndex(field,axis,lower);if(evaluated[face])continue;
            const auto value=stressFace(bound,axis,lower);const double v[3]={value.x,value.y,value.z};
            for(int c=0;c<3;++c)traction[n][(size_t)(axis*3+c)*total+lower]=v[c];evaluated[face]=true;
        }
    }
    for(int cell:cells) {
        int i,j,k;field.getIJK(cell,i,j,k);const auto volume=field.Jac(i,j,k);
        if(!std::isfinite(volume) || volume<=0)throw std::runtime_error("Cell equation has invalid inverse volume");
        const Vector3 value(stencil.velocity.read(cell,0),stencil.velocity.read(cell,1),stencil.velocity.read(cell,2));
        Vector3 source;
        for(const auto& term:terms)if(term.point){const auto q=term.point(field,i,j,k,value);source=source+(term.weight==1?q:Vector3(term.weight*q.x,term.weight*q.y,term.weight*q.z));}
        const double sources[3]={source.x,source.y,source.z};
        for(int c=0;c<3;++c) {
            double spatial=0;bool first=true;
            for(std::size_t n=0;n<terms.size();++n)if(terms[n].spatial){const auto& term=terms[n];auto bound=stencil;bound.coefficient=term.coefficient;bound.constantCoefficient=term.constantCoefficient;bound.traction=&traction[n];const double next=term.weight*term.spatial(bound,cell,c);if(first){spatial=next;first=false;}else spatial+=next;}
            const auto index=(size_t)c*total+cell;rhs[index]=volume*spatial+sources[c];
            if(withoutGradient)(*withoutGradient)[index]=rhs[index];
        }
        for(const auto& term:gradients) {
            const auto g=term.kernel(field,*term.input,cell);
            for(int c=0;c<3;++c)rhs[(size_t)c*total+cell]+=term.weight/term.divisor*g[c];
        }
        for(int c=0;c<3;++c)if(!std::isfinite(rhs[(size_t)c*total+cell]))throw std::runtime_error("Nonfinite cell RHS");
    }
}

inline std::vector<double> candidateFlux(const State::DistributedFieldView& value,const Field& geometry,const std::vector<int>& cells,const std::array<bool,3>& active) {
    std::vector<double> flux((size_t)geometry.TotalSize()*3,0.0);
    for(int cell:cells)for(int axis=0;axis<3;++axis)if(active[axis])for(int sign:{-1,1}) {
        const auto lower=sign>0?cell:detail::adjacent(geometry,cell,axis,-1);
        const auto upper=detail::adjacent(geometry,lower,axis,1);const auto area=detail::faceArea(geometry,axis,lower);
        double q=0;for(int c=0;c<3;++c)q+=0.5*(value.read(lower,c)+value.read(upper,c))*area[c];
        flux[detail::faceIndex(geometry,axis,lower)]=q;
    }return flux;
}

} // namespace SF::Discretization::CellEquation
