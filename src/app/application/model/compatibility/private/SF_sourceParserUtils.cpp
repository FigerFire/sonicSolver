/// @file SF_sourceParserUtils.cpp
/// @brief case 字典解析器的内部无求解逻辑辅助实现。

#include "private/SF_sourceParserUtils.h"

#include "SF_phaseChange.h"

#include <algorithm>
#include <cctype>

namespace SF::IOPrivate {

std::string sourceTokenKey(std::string value) {
    value.erase(std::remove_if(value.begin(), value.end(),
                               [](unsigned char c) {
                                   return c == '_' || c == '-' || c == '+'
                                       || c == ',' || c == ';'
                                       || c == '(' || c == ')'
                                       || std::isspace(c);
                               }),
                value.end());
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return (char)std::tolower(c); });
    return value;
}

bool sourceSchemeMentions(const std::string& scheme,
                          const std::string& token) {
    const std::string key = sourceTokenKey(scheme);
    const std::string needle = sourceTokenKey(token);
    return !needle.empty() && key.find(needle) != std::string::npos;
}

void appendSourceSchemeToken(std::string& scheme, const std::string& token) {
    if (sourceSchemeMentions(scheme, token)) return;
    const std::string current = sourceTokenKey(scheme);
    if (current.empty() || current == "none" || current == "off"
        || current == "false") {
        scheme = token;
        return;
    }
    scheme += "+" + token;
}

bool rpiPhaseChangeActive(
        const Physics::Multiphase::MultiPhaseConfig& config) {
    const std::string model =
        Physics::PhaseChange::normalizeModel(config.phaseChange.model);
    return config.phaseChange.enabled
        && (model == "rpi" || model == "wallboiling");
}

} // namespace SF::IOPrivate
