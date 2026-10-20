#pragma once

/// @file SF_sourceParserUtils.h
/// @brief Source and phase-change dictionary helper declarations.

#include "SF_multiphase.h"

#include <string>

namespace SF::IOPrivate {

std::string sourceTokenKey(std::string value);
bool sourceSchemeMentions(const std::string& scheme,
                          const std::string& token);
/// @brief 向 case-local source 选择串追加一个 token；串由调用方持有。
///
/// scratch 由 CaseAdapter 持有，因此解码过程不依赖任何 parser 全局变量。
void appendSourceSchemeToken(std::string& scheme, const std::string& token);
bool rpiPhaseChangeActive(
    const Physics::Multiphase::MultiPhaseConfig& config);

} // namespace SF::IOPrivate
