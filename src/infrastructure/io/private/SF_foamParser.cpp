/// @file SF_foamParser.cpp
/// @brief case 字典解析器的内部无求解逻辑辅助实现。

#include "private/SF_foamParser.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>

namespace SF::IOPrivate {

std::string foamTrim(std::string s) {
    auto notSpace = [](unsigned char c) { return !std::isspace(c); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), notSpace));
    s.erase(std::find_if(s.rbegin(), s.rend(), notSpace).base(), s.end());
    return s;
}

std::string foamUnquote(std::string s) {
    s = foamTrim(std::move(s));
    while (!s.empty() && (s.back() == ';' || s.back() == ',')) s.pop_back();
    s = foamTrim(std::move(s));
    if (s.size() >= 2 &&
        ((s.front() == '"' && s.back() == '"') ||
         (s.front() == '\'' && s.back() == '\''))) {
        s = s.substr(1, s.size() - 2);
    }
    return foamTrim(std::move(s));
}

std::string foamLower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return (char)std::tolower(c); });
    return s;
}

std::string readTextFile(const std::string& path) {
    std::ifstream in(path);
    if (!in) return {};
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

std::string stripFoamComments(const std::string& text) {
    std::string out;
    out.reserve(text.size());
    bool inBlock = false;
    for (size_t i = 0; i < text.size(); ++i) {
        if (inBlock) {
            if (text[i] == '*' && i + 1 < text.size() && text[i + 1] == '/') {
                inBlock = false;
                ++i;
            }
            continue;
        }
        if (text[i] == '/' && i + 1 < text.size()) {
            if (text[i + 1] == '/') {
                while (i < text.size() && text[i] != '\n') ++i;
                if (i < text.size()) out.push_back('\n');
                continue;
            }
            if (text[i + 1] == '*') {
                inBlock = true;
                ++i;
                continue;
            }
        }
        out.push_back(text[i]);
    }
    return out;
}

bool isFoamWordChar(char c) {
    return std::isalnum((unsigned char)c) || c == '_' || c == '-' || c == '.';
}

size_t skipFoamSpace(const std::string& text, size_t pos) {
    while (pos < text.size() &&
           std::isspace((unsigned char)text[pos])) {
        ++pos;
    }
    return pos;
}

size_t findFoamWord(const std::string& text,
                    const std::string& key,
                    size_t start = 0) {
    for(size_t pos=start;pos<text.size();++pos) {
        if(text[pos]=='"'||text[pos]=='\'') {
            const char quote=text[pos++];
            while(pos<text.size()&&text[pos]!=quote) {
                if(text[pos]=='\\'&&pos+1<text.size())++pos;
                ++pos;
            }
            continue;
        }
        if(text.compare(pos,key.size(),key)!=0)continue;
        const bool beforeOk=pos==0||!isFoamWordChar(text[pos-1]);
        const size_t after=pos+key.size();
        const bool afterOk=after>=text.size()||!isFoamWordChar(text[after]);
        if(beforeOk&&afterOk)return pos;
    }
    return std::string::npos;
}

std::string readFoamWord(const std::string& text, size_t* pos) {
    size_t p = skipFoamSpace(text, *pos);
    if (p >= text.size()) {
        *pos = p;
        return {};
    }
    if (text[p] == '"' || text[p] == '\'') {
        const char quote = text[p++];
        const size_t begin = p;
        while (p < text.size() && text[p] != quote) ++p;
        std::string word = text.substr(begin, p - begin);
        *pos = (p < text.size()) ? p + 1 : p;
        return word;
    }
    const size_t begin = p;
    while (p < text.size() && isFoamWordChar(text[p])) ++p;
    *pos = p;
    return text.substr(begin, p - begin);
}

size_t findMatchingBrace(const std::string& text, size_t open) {
    if (open >= text.size() || text[open] != '{') return std::string::npos;
    int depth = 0;
    for (size_t i = open; i < text.size(); ++i) {
        if (text[i] == '{') ++depth;
        else if (text[i] == '}') {
            --depth;
            if (depth == 0) return i;
        }
    }
    return std::string::npos;
}

std::string foamBlockAfterKey(const std::string& text,
                              const std::string& key) {
    size_t searchFrom = 0;
    while (true) {
        const size_t keyPos = findFoamWord(text, key, searchFrom);
        if (keyPos == std::string::npos) return {};

        size_t p = skipFoamSpace(text, keyPos + key.size());
        if (p < text.size() && text[p] == '{') {
            const size_t close = findMatchingBrace(text, p);
            if (close != std::string::npos) {
                return text.substr(p + 1, close - p - 1);
            }
        }

        searchFrom = keyPos + key.size();
    }
}

std::string foamValueAfterKey(const std::string& text,
                              const std::string& key) {
    const size_t keyPos = findFoamWord(text, key);
    if (keyPos == std::string::npos) return {};
    size_t p = skipFoamSpace(text, keyPos + key.size());
    if (p >= text.size() || text[p] == '{') return {};

    int paren = 0;
    int bracket = 0;
    const size_t begin = p;
    for (; p < text.size(); ++p) {
        if (text[p] == '(') ++paren;
        else if (text[p] == ')' && paren > 0) --paren;
        else if (text[p] == '[') ++bracket;
        else if (text[p] == ']' && bracket > 0) --bracket;
        else if (text[p] == ';' && paren == 0 && bracket == 0) {
            return foamTrim(text.substr(begin, p - begin));
        }
    }
    return foamTrim(text.substr(begin));
}

std::vector<std::string> foamValuesAfterKey(const std::string& text,
                                            const std::string& key) {
    std::vector<std::string> values;
    size_t searchFrom = 0;
    while (true) {
        const size_t keyPos = findFoamWord(text, key, searchFrom);
        if (keyPos == std::string::npos) return values;
        size_t p = skipFoamSpace(text, keyPos + key.size());
        if (p >= text.size() || text[p] == '{') {
            searchFrom = keyPos + key.size();
            continue;
        }

        int paren = 0;
        int bracket = 0;
        const size_t begin = p;
        for (; p < text.size(); ++p) {
            if (text[p] == '(') ++paren;
            else if (text[p] == ')' && paren > 0) --paren;
            else if (text[p] == '[') ++bracket;
            else if (text[p] == ']' && bracket > 0) --bracket;
            else if (text[p] == ';' && paren == 0 && bracket == 0) {
                values.push_back(foamTrim(text.substr(begin, p - begin)));
                break;
            }
        }
        if (p >= text.size()) {
            values.push_back(foamTrim(text.substr(begin)));
            return values;
        }
        searchFrom = p + 1;
    }
}

std::vector<std::pair<std::string, std::string>>
foamChildBlocks(const std::string& block) {
    std::vector<std::pair<std::string, std::string>> children;
    size_t p = 0;
    while (p < block.size()) {
        const std::string name = readFoamWord(block, &p);
        if (name.empty()) {
            ++p;
            continue;
        }
        p = skipFoamSpace(block, p);
        if (p >= block.size() || block[p] != '{') {
            while (p < block.size() && block[p] != ';' && block[p] != '\n') ++p;
            if (p < block.size()) ++p;
            continue;
        }
        const size_t close = findMatchingBrace(block, p);
        if (close == std::string::npos) break;
        children.push_back({name, block.substr(p + 1, close - p - 1)});
        p = close + 1;
    }
    return children;
}

std::vector<double> foamNumbers(std::string text) {
    text = foamTrim(std::move(text));
    if (foamLower(text).rfind("uniform", 0) == 0) {
        text = foamTrim(text.substr(7));
    }
    for (char& c : text) {
        if (c == '(' || c == ')' || c == '[' || c == ']' || c == ',') {
            c = ' ';
        }
    }
    std::stringstream ss(text);
    std::vector<double> values;
    double value = 0.0;
    while (ss >> value) values.push_back(value);
    return values;
}

bool foamScalarValue(const std::string& text, double* value) {
    const std::vector<double> values = foamNumbers(text);
    if (values.empty()) return false;
    *value = values.front();
    return true;
}

bool foamVectorValue(const std::string& text, Vector3* value) {
    const std::vector<double> values = foamNumbers(text);
    if (values.size() < 3) return false;
    *value = Vector3(values[0], values[1], values[2]);
    return true;
}

bool foamBoolValue(const std::string& text, bool fallback) {
    const std::string value = foamLower(foamUnquote(text));
    if (value == "true" || value == "yes" || value == "on" || value == "1") {
        return true;
    }
    if (value == "false" || value == "no" || value == "off" || value == "0") {
        return false;
    }
    return fallback;
}

int foamIntValue(const std::string& text, int fallback) {
    try {
        return std::stoi(foamUnquote(text));
    } catch (...) {
        return fallback;
    }
}

double foamDoubleValue(const std::string& text, double fallback) {
    double value = fallback;
    return foamScalarValue(text, &value) ? value : fallback;
}

std::vector<std::string> foamWordList(std::string text) {
    for (char& c : text) {
        if (c == '(' || c == ')' || c == '[' || c == ']' ||
            c == ',' || c == ';') {
            c = ' ';
        }
    }
    std::stringstream ss(text);
    std::vector<std::string> values;
    std::string token;
    while (ss >> token) {
        token = foamUnquote(token);
        if (!token.empty()) values.push_back(token);
    }
    return values;
}

void appendUnique(std::vector<std::string>& values, const std::string& value) {
    if (value.empty()) return;
    if (std::find(values.begin(), values.end(), value) == values.end()) {
        values.push_back(value);
    }
}

bool foamIntegerTriple(const std::string& text,
                       std::array<int, 3>* values) {
    const std::vector<std::string> tokens = foamWordList(text);
    if (tokens.size() != 3) return false;
    try {
        *values = {{std::stoi(tokens[0]),
                    std::stoi(tokens[1]),
                    std::stoi(tokens[2])}};
        return (*values)[0] > 0 && (*values)[1] > 0 && (*values)[2] > 0;
    } catch (...) {
        return false;
    }
}

BCType foamBCTypeValue(const std::string& text, BCType fallback) {
    const std::string type = foamLower(foamUnquote(text));
    if (type == "fixedvalue" || type == "fixed_value" || type == "noslip") {
        return FIXED_VALUE;
    }
    if (type == "zerogradient" || type == "zero_gradient") {
        return ZERO_GRADIENT;
    }
    if (type == "empty") return EMPTY;
    if (type == "symmetry" || type == "symmetryplane") return SYMMETRY;
    return fallback;
}

ThermalBCType foamThermalBCTypeValue(const std::string& text,
                                     ThermalBCType fallback) {
    std::string type = foamLower(foamUnquote(text));
    type.erase(std::remove_if(type.begin(), type.end(),
                              [](unsigned char c) {
                                  return c == '_' || c == '-' || std::isspace(c);
                              }),
               type.end());
    if (type == "fixedvalue" || type == "fixedtemperature") {
        return ThermalBCType::FixedTemperature;
    }
    if (type == "zerogradient") return ThermalBCType::ZeroGradient;
    if (type == "adiabatic" || type == "insulated") {
        return ThermalBCType::Adiabatic;
    }
    if (type == "heatflux" || type == "fixedheatflux"
        || type == "externalwallheatflux") {
        return ThermalBCType::HeatFlux;
    }
    if (type == "empty") return ThermalBCType::Empty;
    return fallback;
}

bool foamWallHeatSourceType(const std::string& text) {
    std::string type = foamLower(foamUnquote(text));
    type.erase(std::remove_if(type.begin(), type.end(),
                              [](unsigned char c) {
                                  return c == '_' || c == '-' || std::isspace(c);
                              }),
               type.end());
    return type == "wallheat" || type == "wallheatflux"
        || type == "wallheatsource" || type == "heatsource"
        || type == "fixedwallheat";
}

} // namespace SF::IOPrivate
