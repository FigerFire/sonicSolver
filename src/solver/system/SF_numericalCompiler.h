#pragma once

/// @file SF_numericalCompiler.h
/// @brief 将 executable mathematical terms 绑定到 built-in numerical recipes。

#include "SF_numericalSystem.h"
#include "SF_termProviderCatalog.h"
#include "SF_methodObjects.h"

#include <vector>

namespace SF::System::NumericalCompiler {

/// @brief 输入归一化时解析一次，不能从 RHS 重新读 SolverConfig。
void resolveConservativeInputs(const FDM::SolverConfig&,CompiledNumericalSystem&);
const ConservativeExecutionInputs& requireConservativeInputs(const CompiledNumericalSystem&);
CompiledNumericalSystem compile(
    const ExecutableEquationSystem& equations,
    const CompiledExecutionProgram& program,
    const FDM::NumericalRecipeSet& recipes,
    const TermProviderCatalog& catalog,
    const CompiledTimeRecipe& time,
    const ProviderRegistry& providers,
    const std::vector<FormulaOperatorBinding>& bindings = {});

const FDM::TermRecipe& requireUniqueRecipe(
    const CompiledNumericalSystem& system, FDM::TermRole role);
const FDM::TermRecipe* findUniqueRecipe(
    const CompiledNumericalSystem& system, FDM::TermRole role);

} // namespace SF::System::NumericalCompiler
