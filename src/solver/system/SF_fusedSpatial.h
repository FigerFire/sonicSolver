#pragma once
/// @brief 既有融合流体/压力后端的 AST 能力与空间 lowering；不定义公共数学规则。
#include "SF_numericalSystem.h"
#include "SF_termProviderCatalog.h"
namespace SF::System::FusedSpatial {
void lower(const ExecutableEquationSystem&,const CompiledEquationCall&,
    const FDM::NumericalRecipeSet&,const TermProviderCatalog&,
    const std::vector<FormulaOperatorBinding>&,CompiledNumericalSystem&);
void registerProviders(TermProviderCatalog&);
}
