#pragma once

/// @file SF_numericalCompiler.h
/// @brief 将 executable mathematical terms 绑定到 built-in numerical recipes。

#include "SF_resolvedSimulationSystem.h"
#include "SF_termProviderCatalog.h"
#include "SF_methodObjects.h"
#include "SF_numericalSelection.h"
#include "core/system/SF_planFragment.h"

#include <vector>

namespace SF::System { struct BuildRequest; }

namespace SF::System::NumericalCompiler {

/// @brief Resolve typed targets, bind methods/operators and lower the ordered program.
void compileSystem(ResolvedSimulationSystem& result,const ExecutionProgram& program,
    const NumericalSelection& numerics,const std::vector<LegacyPlanFragment>& legacyFragments,
    bool additionalContributions = false);

CompiledNumericalSystem compile(
    const ExecutableEquationSystem& equations,
    const CompiledExecutionProgram& program,
    const FDM::NumericalRecipeSet& recipes,
    const TermProviderCatalog& catalog,
    const CompiledTimeRecipe& time,
    const std::vector<FormulaOperatorBinding>& bindings = {});

const FDM::TermRecipe& requireUniqueRecipe(
    const CompiledNumericalSystem& system, FDM::TermRole role);
const FDM::TermRecipe* findUniqueRecipe(
    const CompiledNumericalSystem& system, FDM::TermRole role);

std::vector<FDM::SourceKind> sourceKinds(
    const CompiledNumericalSystem& system);

} // namespace SF::System::NumericalCompiler
