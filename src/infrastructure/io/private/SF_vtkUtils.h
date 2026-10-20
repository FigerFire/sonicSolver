#pragma once

/// @file SF_vtkUtils.h
/// @brief VTK writer translation units shared private helpers.

#include "SF_resultWriter.h"

#include <algorithm>

namespace SF::IOPrivate {
inline std::string xmlEscape(const std::string& input) {
    std::string output;
    for(char c:input) {
        if(c=='&')output+="&amp;";else if(c=='<')output+="&lt;";
        else if(c=='>')output+="&gt;";else if(c=='"')output+="&quot;";
        else output+=c;
    }
    return output;
}


inline std::string fileNameOnly(const std::string& path) {
    const std::size_t pos = path.find_last_of("/\\");
    return pos == std::string::npos ? path : path.substr(pos + 1);
}

inline std::string outputPrefix(
        const std::string& jobName,
        const std::string& caseName) {
    return jobName.empty() ? caseName : jobName;
}

inline bool hasTemperatureBoundaryScalar(
        const std::vector<ResultWriter::ScalarField>& extraScalars) {
    return std::any_of(
        extraScalars.begin(), extraScalars.end(),
        [](const ResultWriter::ScalarField& scalar) {
            return scalar.writesTemperatureBoundaryId;
        });
}

template <typename T>
void appendBoundaryZoneNames(
        const std::vector<BCSetting<T>>& settings,
        std::vector<std::string>& names) {
    for (const auto& bc : settings) {
        if (std::find(names.begin(), names.end(), bc.name) == names.end()) {
            names.push_back(bc.name);
        }
    }
}

inline void appendBoundaryZoneNames(
        const std::vector<ThermalBCSetting>& settings,
        std::vector<std::string>& names) {
    for (const auto& bc : settings) {
        if (std::find(names.begin(), names.end(), bc.name) == names.end()) {
            names.push_back(bc.name);
        }
    }
}

} // namespace SF::IOPrivate
