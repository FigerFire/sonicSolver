#pragma once
#include "SF_methodObjects.h"
#include "solver/discretization/SF_cellEquation.h"

namespace SF::System {
/// @brief 已完整消费的 cell 向量调用响应；索引引用 operators，不保存第二份数学。
struct CompiledCellGradient {
    std::string input;double weight=1,divisor=1;
    Discretization::CellEquation::GradientKernel kernel=nullptr;
};
struct CompiledCellEquation {
    std::string potential,capacity;
    double capacityValue=1.0,diffusivity=0.0,advectionWeight=0.0;
    std::vector<std::size_t> terms;
    std::vector<CompiledCellGradient> gradients;
};

inline const CompiledCellEquation* cellEquation(const CompiledEquationCall& call) {
    return std::any_cast<CompiledCellEquation>(&call.providerContract);
}

CompiledEquationCall compileCellEquation(const ExecutableEquationSystem&,const EquationCall&,const NumericalBinding&);
void lowerCellEquation(const ExecutableEquationSystem&,CompiledEquationCall&,const FDM::NumericalRecipeSet&,
    const TermProviderCatalog&,const std::vector<FormulaOperatorBinding>&,CompiledNumericalSystem&);
}
