#pragma once
/// @file SF_eulerianRelations.h
/// @brief 现有双欧拉内核的数学关系；不含时间步或耦合顺序。
#include "core/system/SF_equationIR.h"
namespace SF::System {
std::vector<Equation> eulerianPhaseRelations(const std::string& phase,bool reference);
std::vector<Equation> eulerianPressureRelations();
/// @brief 仅投影旧 assembler 检查的 term inventory；AST 保留完整符号和符号运算。
SF::Equation::Definition eulerianBackendDefinition(const Equation& equation);
void projectEulerianBackendDefinitions(ExecutableEquationSystem& system);
}
