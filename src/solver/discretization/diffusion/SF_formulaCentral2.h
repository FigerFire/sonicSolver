#pragma once

/// @file SF_formulaCentral2.h
/// @brief Cartesian Central2 provider for div(nu * grad(X)).

#include "solver/discretization/SF_formulaOperator.h"

#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

namespace SF::Discretization {

struct FormulaCartesianGrid {
    int nx = 0;
    int ny = 0;
    double dx = 0.0;
    double dy = 0.0;
};

namespace detail {

inline bool central2Pattern(const System::FormulaExpr& expression) {
    using Kind=System::FormulaExpr::Kind;
    if (expression.kind!=Kind::Operator || expression.name!="div"
        || expression.arguments.size()!=1) return false;
    const auto& product=expression.arguments[0];
    if (product.kind!=Kind::Multiply || product.arguments.size()!=2
        || product.arguments[0].kind!=Kind::Symbol) return false;
    const auto& gradient=product.arguments[1];
    return gradient.kind==Kind::Operator && gradient.name=="grad"
        && gradient.arguments.size()==1
        && gradient.arguments[0].kind==Kind::Symbol;
}

inline std::pair<std::string,std::string> central2Symbols(
        const System::FormulaExpr& expression) {
    if (!central2Pattern(expression))
        throw std::runtime_error("Central2 diffusion requires div(nu * grad(X)).");
    const auto& product=expression.arguments[0];
    return {product.arguments[0].name,product.arguments[1].arguments[0].name};
}

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
inline System::FormulaOperatorProvider central2FormulaDiffusion(
        FormulaCartesianGrid grid) {
    System::FormulaOperatorProvider provider;
    provider.id="central2.cartesian.diffusion";
    provider.mathematicalOperator="div";
    provider.matches=[](const System::FormulaExpr& expression) {
        return detail::central2Pattern(expression);
    };
    provider.compileValue=[grid](const System::FormulaExpr& expression) {
        const auto symbols=detail::central2Symbols(expression);
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
    };
    provider.compileLinear=[grid](const System::FormulaExpr& expression,
                                   const std::string& target) {
        const auto symbols=detail::central2Symbols(expression);
        if (symbols.second!=target)
            throw std::runtime_error("Central2 implicit diffusion target does not "
                                     "match grad(X).");
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
    };
    return provider;
}

} // namespace SF::Discretization
