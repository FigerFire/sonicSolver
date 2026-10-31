#pragma once
/// @file SF_declarations.h
/// @brief 普通 native 输入解码为现役四轴 IR，不解释数值算法。
#include "core/model/SF_model.h"
#include "core/system/SF_systemContribution.h"
namespace SF::System::Declarations {
/// @brief 与 native WHAT 共用的 typed AST 解码；symbols 必须已声明。
FormulaExpr formula(const Model::Parameters&,const std::string&,const FormulaSymbols&,int depth=0);
/// @brief 普通 native HOW 解码；内部算法直接构造 typed IR，不经过此入口。
ExecutionProgram execution(const Model::Parameters&,int formatVersion,Provenance);
/// @brief 解码方法选择；默认/显式优先级由调用方按来源冻结。
std::vector<NumericalBinding> numerics(const Model::Parameters&,Provenance);
}
