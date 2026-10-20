#pragma once

/// @file SF_valueParser.h
/// @brief 已解析 native 参数的数值/列表转换；不读取文件或解析字典。

#include "core/state/SF_valueTypes.h"

#include <array>
#include <string>
#include <utility>
#include <vector>

namespace SF::IOPrivate {

std::string valueTrim(std::string value);
std::string valueUnquote(std::string value);
std::string valueLower(std::string value);
std::string readValueWord(const std::string& text,std::size_t* position);
std::vector<double> valueNumbers(std::string text);
bool valueScalarValue(const std::string& text, double* value);
bool valueVectorValue(const std::string& text, Vector3* value);
bool valueBoolValue(const std::string& text, bool fallback);
int valueIntValue(const std::string& text, int fallback);
double valueDoubleValue(const std::string& text, double fallback);
std::vector<std::string> valueWordList(std::string text);
void appendUnique(std::vector<std::string>& values,
                  const std::string& value);
bool valueIntegerTriple(const std::string& text,
                       std::array<int, 3>* values);
/// @brief 解析边界约束；slip/noSlip 仅允许在速度边界，展开为已有中立约束。
/// @param velocityBoundary true 表示速度 patch，不用于初始条件或标量。
BCType valueBCTypeValue(const std::string& text, BCType fallback,
                      bool velocityBoundary = false);
ThermalBCType valueThermalBCTypeValue(
    const std::string& text,
    ThermalBCType fallback);
bool valueWallHeatSourceType(const std::string& text);

} // namespace SF::IOPrivate
