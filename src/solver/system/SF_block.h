#pragma once
#include "SF_methodObjects.h"
#include "solver/discretization/SF_formulaOperator.h"
namespace SF::System {
/// @brief 串行标量块的不可变AST lowering；行标签来自唯一主扩散目标。
struct CompiledLinearBlock {
    std::vector<std::string> equations,targets,rowTargets;
    std::vector<StateUse> valueReads;
    std::vector<FormulaLinearKernel> rows;
    FDM::LinearSolverConfig linear;
};
/// @brief 已实现的稳态 Central2 块装配及 HYPRE 方法。
const IProvider& linearBlockMethod();
}
