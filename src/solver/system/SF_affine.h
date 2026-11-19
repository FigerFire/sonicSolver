#pragma once
#include "SF_methodObjects.h"
#include "solver/discretization/SF_formulaOperator.h"
namespace SF::System {
/// @brief 针对HOW目标的局部仿射数值实现；WHAT仍保存原等式。
struct CompiledAffineEquation {
    WorkspaceBoundary boundaryPreparation=WorkspaceBoundary::None;
    FormulaValueKernel coefficient,remainder;
};
/// @brief 以 A*target+R 分解局部残量，冻结结果发布动作。
CompiledEquationCall compileLocalAffine(const ExecutableEquationSystem&,const EquationCall&,const NumericalBinding&);
}
