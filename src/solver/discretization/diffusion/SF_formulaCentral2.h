#pragma once

/// @file SF_formulaCentral2.h
/// @brief Cartesian Central2 provider for div(nu * grad(X)).

#include "solver/discretization/SF_formulaOperator.h"

#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

namespace SF::Discretization {

/// Constant-coefficient nodal Central2. Geometry and stage neighbors are
/// frozen runtime bindings; this operator never selects a temporal method.
/// @brief 借用已发布的输入分量；交叉扩散不复制 partner state。
inline System::FormulaValueKernel central2StageKernel(double D,std::size_t slot,bool scalar) {
    return [D,slot,scalar](const System::FormulaValues& v,int cell,int component) {
        if(scalar)component=0;
        double result=0;
        for(int axis=0;axis<3;++axis) {
            const double h=v.spacing(axis);if(h==0)continue;
            result+=D*(v.boundNeighbor(slot,cell,axis,-1,component)-2*v.boundNeighbor(slot,cell,axis,0,component)
                +v.boundNeighbor(slot,cell,axis,1,component))/(h*h);
        }
        return result;
    };
}
inline System::FormulaValueKernel conservativeCentral2Kernel(double D,int coefficientSlot,bool primitive) {
    return [D,coefficientSlot,primitive](const auto& v,int cell,int) {
            double result=0;
            for(int axis=0;axis<3;++axis) {
                const double h=v.spacing(axis);if(h==0)continue;
                const double rho=v.boundReads[1](cell,0);
                if(!std::isfinite(rho) || rho<=0)throw std::runtime_error("Scalar transport Stage density is not positive.");
                const double center=v.boundReads[0](cell,0)/(primitive?1:rho);
                for(int sign:{-1,1}) {
                    const double rn=v.boundNeighbor(1,cell,axis,sign,0);
                    if(!std::isfinite(rn) || rn<=0)throw std::runtime_error("Scalar transport neighbor density is not positive.");
                    const double cn=v.boundNeighbor(0,cell,axis,sign,0)/(primitive?1:rn);
                    const double centerD=coefficientSlot<0?D:v.boundReads.at(coefficientSlot)(cell,0);
                    const double nextD=coefficientSlot<0?D:v.boundNeighbor(coefficientSlot,cell,axis,sign,0);
                    if(!std::isfinite(centerD) || !std::isfinite(nextD) || centerD<0 || nextD<0)
                        throw std::runtime_error("Scalar Stage diffusion coefficient must be finite and nonnegative.");
                    result+=(coefficientSlot<0?.5*D*(rho+rn):.5*(centerD*rho+nextD*rn))*(cn-center)/(h*h);
                }
            }
            return result;
        };
}
/// A single face expression is reused with +/- divergence signs. This is an
/// explicit scalar Upwind1 recipe, not a replacement for the NS flux recipe.
inline double conservativeScalarAdvection(const System::FormulaValues& v,int cell) {
    double divergence=0;
    for(int axis=0;axis<3;++axis) {
        const double h=v.spacing(axis);if(h==0)continue;
        const double rho=v.boundReads[1](cell,0),m=v.boundReads[2](cell,axis);
        if(!std::isfinite(rho) || rho<=0)throw std::runtime_error("Scalar transport density is not positive.");
        const double c=v.boundReads[0](cell,0)/rho;
        for(int sign:{-1,1}) {
            const double rn=v.boundNeighbor(1,cell,axis,sign,0),mn=v.boundNeighbor(2,cell,axis,sign,axis);
            if(!std::isfinite(rn) || rn<=0)throw std::runtime_error("Scalar transport neighbor density is not positive.");
            const double cn=v.boundNeighbor(0,cell,axis,sign,0)/rn;
            const double mass=.5*(m+mn);
            const double left=sign>0?c:cn,right=sign>0?cn:c;
            divergence+=sign*mass*(mass>=0?left:right)/h;
        }
    }return divergence;
}

struct FormulaCartesianGrid {
    int nx = 0;
    int ny = 0;
    double dx = 0.0;
    double dy = 0.0;
};

namespace detail {

template<class Contribution>
void central2Faces(const FormulaCartesianGrid& grid,
                   const System::FormulaValues& values,
                   const std::string& coefficient,
                   const std::string& target,
                   int cell,int component,Contribution&& contribute) {
    if (grid.nx<1 || grid.ny<1 || !(grid.dx>0.0)
        || !(grid.dy>0.0) || cell<0 || cell>=grid.nx*grid.ny
        || values.ownedCells!=grid.nx*grid.ny || !values.boundaryValue)
        throw std::runtime_error("Central2 Equation has invalid Cartesian geometry.");
    const int x=cell%grid.nx;
    const int y=cell/grid.nx;
    const double center=values.read(coefficient,cell,0);
    if (!std::isfinite(center) || !(center>0.0))
        throw std::runtime_error("Central2 diffusion coefficient is not positive and finite.");
    for (int axis=0;axis<2;++axis) {
        if (axis==1 && grid.ny==1) continue;
        const double spacing=axis==0?grid.dx:grid.dy;
        for (int sign:{-1,1}) {
            const int nextX=x+(axis==0?sign:0);
            const int nextY=y+(axis==1?sign:0);
            const bool interior=nextX>=0 && nextX<grid.nx
                && nextY>=0 && nextY<grid.ny;
            const int next=interior?nextY*grid.nx+nextX:-1;
            const double neighbourNu=interior
                ? values.read(coefficient,next,0) : center;
            if (!std::isfinite(neighbourNu) || !(neighbourNu>0.0))
                throw std::runtime_error("Central2 face coefficient is invalid.");
            const double faceNu=0.5*(center+neighbourNu);
            contribute(axis,sign,next,faceNu/(spacing*spacing),target);
        }
    }
}

} // namespace detail

/// @brief Bind a grid and the existing state views once; the compiled provider
/// contributes coefficients or values without choosing a solver lifecycle.
inline System::FormulaValueKernel central2CartesianValue(FormulaCartesianGrid grid,std::pair<std::string,std::string> symbols) {
    return System::FormulaValueKernel{
            [grid,symbols](const System::FormulaValues& values,int cell,int component) {
                const double center=values.read(symbols.second,cell,component);
                double sum=0.0;
                detail::central2Faces(grid,values,symbols.first,symbols.second,
                    cell,component,[&](int axis,int sign,int next,double coefficient,
                                        const std::string& target) {
                        const double neighbour=next>=0
                            ? values.read(target,next,component)
                            : values.boundaryValue(target,cell,axis,sign,component);
                        sum+=coefficient*(neighbour-center);
                    });
                return sum;
            }};
}
inline System::FormulaLinearKernel central2CartesianLinear(FormulaCartesianGrid grid,std::pair<std::string,std::string> symbols) {
    return System::FormulaLinearKernel{
            [grid,symbols](const System::FormulaValues& values,int cell,int component) {
                System::FormulaLinearTerm term;
                detail::central2Faces(grid,values,symbols.first,symbols.second,
                    cell,component,[&](int axis,int sign,int next,double coefficient,
                                        const std::string& target) {
                        term.add(values.dof(target,cell,component),-coefficient);
                        if (next>=0)
                            term.add(values.dof(target,next,component),coefficient);
                        else
                            term.constant+=coefficient*values.boundaryValue(
                                target,cell,axis,sign,component);
                    });
                return term;
            }};
}

} // namespace SF::Discretization
