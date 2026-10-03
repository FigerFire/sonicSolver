#pragma once

/// @file SF_formulaCompiler.h
/// @brief Compile HOW FormulaCalls using numerical-provider capabilities.

#include "core/system/SF_formula.h"
#include "core/system/SF_solveProgram.h"
#include "solver/discretization/SF_formulaOperator.h"
#include "SF_termProviderCatalog.h"

#include <functional>
#include <string>
#include <vector>

namespace SF::System {


namespace Legacy {
/// @brief Isolated legacy FormulaCompiler request mode. Production methods
/// derive realization from their concrete method identity.
enum class FormulaMode { Assign, Explicit, Implicit };

struct FormulaCall {
    std::string equation;
    std::string target;
    FormulaMode mode = FormulaMode::Assign;
};

} // namespace Legacy

class FormulaOperatorCatalog {
public:
    void add(FormulaOperatorProvider provider);
    const FormulaOperatorProvider& resolve(
        const Equation& formula,const FormulaExpr& expression,
        const std::vector<FormulaOperatorBinding>& bindings) const;
private:
    std::vector<FormulaOperatorProvider> providers_;
};

/// @brief Frozen execution binding. Solve selection is supplied by the caller;
/// this object only produces a GlobalDofSystem and publishes its solution.
struct CompiledFormulaCall {
    Legacy::FormulaCall call;
    std::string backend;
    std::vector<std::string> numericalProviders;
    std::function<void(FormulaValues&)> assign;
    /// @brief Explicit PDE derivative for the selected target. Stage/update
    /// coefficients remain owned by CompiledTimeRecipe.
    FormulaValueKernel evaluateRhs;
    std::function<LinearAlgebra::GlobalDofSystem(const FormulaValues&)> assemble;
    void writeSolution(FormulaValues& values,
                       const LinearAlgebra::SolveResult& result) const;
};

class FormulaCompiler {
public:
    static CompiledFormulaCall compile(
        const EquationRegistry& formulas,const Legacy::FormulaCall& call,
        const FormulaOperatorCatalog& providers,
        const std::vector<FormulaOperatorBinding>& bindings);
};

} // namespace SF::System
