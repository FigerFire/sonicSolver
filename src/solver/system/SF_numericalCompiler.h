#pragma once

/// @file SF_numericalCompiler.h
/// @brief 将 executable mathematical terms 绑定到 built-in numerical recipes。

#include "SF_resolvedSimulationSystem.h"

#include <vector>

namespace SF::System::NumericalCompiler {

CompiledNumericalSystem compile(
    const ExecutableEquationSystem& equations,
    const FDM::NumericalRecipeSet& recipes);

const FDM::TermRecipe& requireUniqueRecipe(
    const CompiledNumericalSystem& system, FDM::TermRole role);
const FDM::TermRecipe* findUniqueRecipe(
    const CompiledNumericalSystem& system, FDM::TermRole role);

std::vector<FDM::SourceKind> sourceKinds(
    const CompiledNumericalSystem& system);

} // namespace SF::System::NumericalCompiler
