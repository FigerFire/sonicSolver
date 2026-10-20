#include "private/SF_valueParser.h"
#include <algorithm>
#include <cctype>
#include <sstream>
#include <stdexcept>
namespace SF::IOPrivate {
std::string valueTrim(std::string s) {
    auto notSpace = [](unsigned char c) { return !std::isspace(c); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), notSpace));
    s.erase(std::find_if(s.rbegin(), s.rend(), notSpace).base(), s.end());
    return s;
}

std::string valueUnquote(std::string s) {
    s = valueTrim(std::move(s));
    while (!s.empty() && (s.back() == ';' || s.back() == ',')) s.pop_back();
    s = valueTrim(std::move(s));
    if (s.size() >= 2 &&
        ((s.front() == '"' && s.back() == '"') ||
         (s.front() == '\'' && s.back() == '\''))) {
        s = s.substr(1, s.size() - 2);
    }
    return valueTrim(std::move(s));
}

std::string valueLower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return (char)std::tolower(c); });
    return s;
}


std::string readValueWord(const std::string& text, size_t* pos) {
    size_t p=*pos;while(p<text.size() && std::isspace(static_cast<unsigned char>(text[p])))++p;
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
    while (p < text.size() && (std::isalnum(static_cast<unsigned char>(text[p])) || text[p]=='_' || text[p]=='-' || text[p]=='.' || text[p]=='/')) ++p;
    *pos = p;
    return text.substr(begin, p - begin);
}

std::vector<double> valueNumbers(std::string text) {
    text = valueTrim(std::move(text));
    if (valueLower(text).rfind("uniform", 0) == 0) {
        text = valueTrim(text.substr(7));
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

bool valueScalarValue(const std::string& text, double* value) {
    const std::vector<double> values = valueNumbers(text);
    if (values.empty()) return false;
    *value = values.front();
    return true;
}

bool valueVectorValue(const std::string& text, Vector3* value) {
    const std::vector<double> values = valueNumbers(text);
    if (values.size() < 3) return false;
    *value = Vector3(values[0], values[1], values[2]);
    return true;
}

bool valueBoolValue(const std::string& text, bool fallback) {
    const std::string value = valueLower(valueUnquote(text));
    if (value == "true" || value == "yes" || value == "on" || value == "1") {
        return true;
    }
    if (value == "false" || value == "no" || value == "off" || value == "0") {
        return false;
    }
    return fallback;
}

int valueIntValue(const std::string& text, int fallback) {
    try {
        return std::stoi(valueUnquote(text));
    } catch (...) {
        return fallback;
    }
}

double valueDoubleValue(const std::string& text, double fallback) {
    double value = fallback;
    return valueScalarValue(text, &value) ? value : fallback;
}

std::vector<std::string> valueWordList(std::string text) {
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
        token = valueUnquote(token);
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

bool valueIntegerTriple(const std::string& text,
                       std::array<int, 3>* values) {
    const std::vector<std::string> tokens = valueWordList(text);
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

BCType valueBCTypeValue(const std::string& text, BCType fallback,
                      bool velocityBoundary) {
    const std::string type = valueLower(valueUnquote(text));
    if (type.empty()) return fallback;
    if (type == "noslip" || type == "slip") {
        if (!velocityBoundary)
            throw std::runtime_error("Boundary law '"+type
                +"' is only valid for velocity patches; scalar, thermal and initial conditions are independent.");
        return type == "noslip" ? FIXED_VALUE : SYMMETRY;
    }
    if (type == "fixedvalue" || type == "fixed_value") {
        return FIXED_VALUE;
    }
    if (type == "zerogradient" || type == "zero_gradient") {
        return ZERO_GRADIENT;
    }
    if (type == "empty") return EMPTY;
    if (type == "symmetry" || type == "symmetryplane") return SYMMETRY;
    throw std::runtime_error("Unknown boundary law '"+type+"'; no boundary fallback is applied.");
}

ThermalBCType valueThermalBCTypeValue(const std::string& text,
                                     ThermalBCType fallback) {
    std::string type = valueLower(valueUnquote(text));
    type.erase(std::remove_if(type.begin(), type.end(),
                              [](unsigned char c) {
                                  return c == '_' || c == '-' || std::isspace(c);
                              }),
               type.end());
    if (type == "noslip" || type == "slip")
        throw std::runtime_error("slip/noSlip are velocity laws; choose the thermal boundary independently.");
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

bool valueWallHeatSourceType(const std::string& text) {
    std::string type = valueLower(valueUnquote(text));
    type.erase(std::remove_if(type.begin(), type.end(),
                              [](unsigned char c) {
                                  return c == '_' || c == '-' || std::isspace(c);
                              }),
               type.end());
    return type == "wallheat" || type == "wallheatflux"
        || type == "wallheatsource" || type == "heatsource"
        || type == "fixedwallheat";
}


}
