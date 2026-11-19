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
/// @brief 编译一个已经定位并解析了 STATE 输入的 HOW 叶节点。保留调用身份及局部效果。
CompiledEquationCall compileCall(const ExecutableEquationSystem&,const StateRegistry&,
    const ExecutionScope&,const std::vector<NumericalBinding>&,
    const ProviderRegistry&,const CompiledTimeRecipe*);

CompiledNumericalSystem compile(
    const ExecutableEquationSystem& equations,
    CompiledExecutionProgram& program,
    const FDM::NumericalRecipeSet& recipes,
    const TermProviderCatalog& catalog,
    const CompiledTimeRecipe& time,
    const std::vector<FormulaOperatorBinding>& bindings = {});

const FDM::TermRecipe& requireUniqueRecipe(
    const CompiledNumericalSystem& system, FDM::TermRole role);
const FDM::TermRecipe* findUniqueRecipe(
    const CompiledNumericalSystem& system, FDM::TermRole role);

} // namespace SF::System::NumericalCompiler
