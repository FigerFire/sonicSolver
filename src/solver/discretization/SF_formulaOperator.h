#pragma once

/// @file SF_formulaOperator.h
/// @brief Numerical-provider contract for Equation AST operator occurrences.

#include "core/system/SF_formula.h"
#include "solver/linearAlgebra/SF_globalDofSystem.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

namespace SF::System {

/// @brief Runtime views supplied by state realization without copying state.
struct FormulaValues {
    /// Compile-time symbol slots, bound once; no per-cell string lookup.
    std::vector<std::function<double(int,int)>> boundReads;
    std::function<double(std::size_t,int,int,int,int)> boundNeighbor;
    /// Frozen orthogonal-stencil geometry/neighbor bindings supplied by a provider.
    std::function<double(int)> spacing;
    std::function<double(int,int,int)> neighbor;
    int ownedCells = 0;
    int components = 1;
    std::function<double(const std::string&,int,int)> read;
    std::function<void(const std::string&,int,int,double)> write;
    std::function<LinearAlgebra::GlobalDofId(const std::string&,int,int)> dof;
    std::function<double(const std::string&,int,int,int,int)> boundaryValue;
};

/// @brief Affine contribution F(X)=sum(a_i X_i)+constant for one row.
struct FormulaLinearTerm {
    std::vector<LinearAlgebra::GlobalDofCoefficient> coefficients;
    double constant = 0.0;
    void add(LinearAlgebra::GlobalDofId column,double value) {
        if (!column.valid() || !std::isfinite(value))
            throw std::runtime_error("Equation linear contribution has invalid coefficient.");
        const auto found=std::find_if(coefficients.begin(),coefficients.end(),
            [&](const auto& item) { return item.column==column; });
        if (found==coefficients.end()) coefficients.push_back({column,value});
        else found->value+=value;
    }
};

using FormulaValueKernel =
    std::function<double(const FormulaValues&,int,int)>;
using FormulaLinearKernel =
    std::function<FormulaLinearTerm(const FormulaValues&,int,int)>;

/// @brief A numerical provider compiles one mathematical operator occurrence.
/// Discretization mathematics lives here, never in FormulaCompiler.
struct FormulaOperatorProvider {
    std::string id;
    std::string mathematicalOperator;
    std::function<bool(const FormulaExpr&)> matches;
    std::function<FormulaValueKernel(const FormulaExpr&)> compileValue;
    std::function<FormulaLinearKernel(const FormulaExpr&,const std::string&)>
        compileLinear;
};

} // namespace SF::System
