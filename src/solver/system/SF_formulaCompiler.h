#pragma once

/// @file SF_formulaCompiler.h
/// @brief Typed structural extraction and numerical value-provider binding.

#include "core/system/SF_formula.h"
#include "core/system/SF_solveProgram.h"
#include "solver/discretization/SF_formulaOperator.h"
#include "SF_termProviderCatalog.h"

#include <functional>
#include <string>
#include <vector>

namespace SF::System {


/// @brief Immutable structural index over the authored AST, not another equation authority.
struct FormulaNode {
    const FormulaExpr* expression=nullptr;
    std::string path;
    unsigned targetDegree=0; // 0: independent; 1: affine candidate; 2: nonlinear/unknown.
};
struct FormulaStructure {
    std::vector<FormulaNode> nodes;
    std::vector<std::string> symbols;
    std::vector<Math::Op> operations;
    unsigned targetDegree=0;
};
FormulaStructure inspectFormula(const FormulaExpr& expression,std::string_view target={});
/// @brief Extract only a single positive additive ddt(target); no mass inversion or term reordering.
FormulaExpr explicitFormulaRhs(const Equation& equation,std::string_view target);

/// @brief A numerical provider compiles one mathematical operator occurrence.
/// Discretization mathematics lives here, never in FormulaCompiler.
struct FormulaOperatorProvider {
    std::string id;
    Math::Op mathematicalOperator=Math::Op::Diffusion;
    std::function<bool(const FormulaExpr&)> matches;
    std::function<FormulaValueKernel(const FormulaExpr&)> compileValue;
    std::function<FormulaLinearKernel(const FormulaExpr&,const std::string&)>
        compileLinear;
    std::string extensionName;
    Math::Signature (*extensionRule)(const std::vector<Math::Signature>&)=nullptr;
};


class FormulaOperatorCatalog {
public:
    void add(FormulaOperatorProvider provider);
    const FormulaOperatorProvider& resolve(
        const Equation& formula,const FormulaExpr& expression,
        const std::vector<FormulaOperatorBinding>& bindings) const;
private:
    std::vector<FormulaOperatorProvider> providers_;
};

/// Shared structural value lowering; no legacy execution request or solve mode.
FormulaValueKernel compileFormulaValue(const Equation& equation,const FormulaExpr& expression,
    const FormulaOperatorCatalog& providers,std::vector<std::string>& used,
    const std::vector<std::string>& boundSymbols={},
    const StateRegistry* state=nullptr,std::vector<std::string>* trace=nullptr);

} // namespace SF::System

namespace SF::Discretization {
struct FormulaCartesianGrid;
// AST binding is compiler-side. Kernels and geometry stay in discretization.
System::FormulaOperatorProvider central2ScalarDiffusion(const System::StateRegistry&,const std::vector<std::string>& slots);
System::FormulaOperatorProvider conservativeScalarCentral2(double D,int coefficientSlot=-1,bool primitive=false);
System::FormulaOperatorProvider central2FormulaDiffusion(FormulaCartesianGrid grid);
}
