#pragma once
/// @file SF_eulerianRelations.h
/// @brief 现有双欧拉内核的数学关系；不含时间步或耦合顺序。
#include "core/system/SF_equationIR.h"
namespace SF::System {
std::vector<Equation> eulerianPhaseRelations(const std::string& phase,bool reference);
std::vector<Equation> eulerianPressureRelations();
}
