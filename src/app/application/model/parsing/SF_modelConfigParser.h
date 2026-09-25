#pragma once

#include "SF_parserCommon.h"
#include "core/config/SF_configTypes.h"
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace SF::FDM {
/// @brief Parse source-term family names.
/// @param value Source list such as `"Gravity+MRF"` or `"None"`.
/// @return Ordered source kinds; disabled/none tokens are omitted.
inline std::vector<SourceKind> parseSourceKinds(const std::string& value) {
    std::vector<SourceKind> kinds;
    std::string sourceText = value.empty() ? "None" : value;
    for (const std::string& raw : splitList(sourceText)) {
        std::string t = normalizeToken(raw);
        if (t.empty() || t == "none" || t == "off" || t == "false") continue;
        if (t == "gravity") kinds.push_back(SourceKind::Gravity);
        else if (t == "mrf" || t == "rotating" || t == "rotation") kinds.push_back(SourceKind::MRF);
        else if (t == "wallheat" || t == "wallheatsource"
                 || t == "wallheatflux" || t == "heat" || t == "heater") {
            kinds.push_back(SourceKind::WallHeat);
        }
        else {
            throw std::invalid_argument(
                "Unsupported source term '" + raw
                + "'. Supported source terms: None, Gravity, MRF, WallHeat.");
        }
    }
    return kinds;
}
/// @brief 解析湍流大类字符串。
/// @param value 用户配置中的 type，例如 `"RAS"`、`"LES"`、`"DNS"`。
/// @return 强类型湍流大类；空值返回 `None`。
inline TurbulenceFamily parseTurbulenceFamily(const std::string& value) {
    std::string t = normalizeToken(value);
    if (t.empty() || t == "none" || t == "off" || t == "false"
        || t == "laminar") {
        return TurbulenceFamily::None;
    }
    if (t == "ras" || t == "rans") return TurbulenceFamily::RAS;
    if (t == "les") return TurbulenceFamily::LES;
    if (t == "dns") return TurbulenceFamily::DNS;
    throw std::invalid_argument(
        "Unsupported turbulence family '" + value
        + "'. Supported families: None, RAS, LES, DNS. Recommended "
          "pairings: RAS+kEpsilon/kOmegaSST, LES+Smagorinsky, DNS+DNS.");
}
/// @brief 解析具体湍流模型字符串。
/// @param value 用户配置中的 value/model，例如 `"kEpsilon"`、`"Smagorinsky"`。
/// @return 强类型模型标识；空值返回 `None`。
inline TurbulenceModelKind parseTurbulenceModelKind(const std::string& value) {
    std::string t = normalizeToken(value);
    if (t.empty() || t == "none" || t == "off" || t == "false"
        || t == "laminar") {
        return TurbulenceModelKind::None;
    }
    if (t == "kepsilon" || t == "kep" || t == "ke") return TurbulenceModelKind::kEpsilon;
    if (t == "komegasst" || t == "sst") return TurbulenceModelKind::kOmegaSST;
    if (t == "smagorinsky" || t == "smago") return TurbulenceModelKind::Smagorinsky;
    if (t == "dns") return TurbulenceModelKind::DNS;
    throw std::invalid_argument(
        "Unsupported turbulence model '" + value
        + "'. Supported models: None, kEpsilon, kOmegaSST, Smagorinsky, DNS. "
          "Recommended pairings: RAS+kEpsilon/kOmegaSST, LES+Smagorinsky, "
          "DNS+DNS.");
}
} // namespace SF::FDM
