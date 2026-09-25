#pragma once

/// @file SF_foamParser.h
/// @brief 旧 case 兼容适配器和 blockMesh/SFM 所需的词法工具。

#include "core/state/SF_valueTypes.h"

#include <array>
#include <string>
#include <utility>
#include <vector>

namespace SF::IOPrivate {

std::string foamTrim(std::string value);
std::string foamUnquote(std::string value);
std::string foamLower(std::string value);
std::string readTextFile(const std::string& path);
std::string stripFoamComments(const std::string& text);
std::size_t skipFoamSpace(const std::string& text, std::size_t pos);
std::string readFoamWord(const std::string& text, std::size_t* pos);
std::size_t findMatchingBrace(
    const std::string& text,
    std::size_t open);
std::string foamBlockAfterKey(const std::string& text,
                              const std::string& key);
std::string foamValueAfterKey(const std::string& text,
                              const std::string& key);
std::vector<std::string> foamValuesAfterKey(
    const std::string& text,
    const std::string& key);
std::vector<std::pair<std::string, std::string>>
foamChildBlocks(const std::string& block);
std::vector<double> foamNumbers(std::string text);
bool foamScalarValue(const std::string& text, double* value);
bool foamVectorValue(const std::string& text, Vector3* value);
bool foamBoolValue(const std::string& text, bool fallback);
int foamIntValue(const std::string& text, int fallback);
double foamDoubleValue(const std::string& text, double fallback);
std::vector<std::string> foamWordList(std::string text);
void appendUnique(std::vector<std::string>& values,
                  const std::string& value);
bool foamIntegerTriple(const std::string& text,
                       std::array<int, 3>* values);
BCType foamBCTypeValue(const std::string& text, BCType fallback);
ThermalBCType foamThermalBCTypeValue(
    const std::string& text,
    ThermalBCType fallback);
bool foamWallHeatSourceType(const std::string& text);

} // namespace SF::IOPrivate
