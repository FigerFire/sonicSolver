#pragma once

#include "core/config/SF_configTypes.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

namespace SF::FDM {
/// @brief Normalize user-facing config tokens for tolerant parsing.
/// @param value Raw token from `.sf`/TOML-like config or legacy caller.
/// @return Lowercase token with quotes, underscores, dashes, and whitespace removed.
inline std::string normalizeToken(std::string value) {
    value.erase(std::remove(value.begin(), value.end(), '"'), value.end());
    value.erase(std::remove_if(value.begin(), value.end(),
                               [](unsigned char c) {
                                   return c == '_' || c == '-' || std::isspace(c);
                               }),
                value.end());
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return (char)std::tolower(c); });
    return value;
}
/// @brief Split a source-list style string.
/// @param value Text using comma, semicolon, plus, or whitespace separators.
/// @return Ordered non-empty tokens; tokens are not normalized here.
inline std::vector<std::string> splitList(const std::string& value) {
    std::vector<std::string> out;
    std::string token;
    for (char c : value) {
        if (c == ',' || c == ';' || c == '+' || std::isspace((unsigned char)c)) {
            if (!token.empty()) {
                out.push_back(token);
                token.clear();
            }
        } else {
            token.push_back(c);
        }
    }
    if (!token.empty()) out.push_back(token);
    return out;
}
} // namespace SF::FDM
