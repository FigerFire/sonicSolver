/// @file SF_sourceContribution.cpp
/// @brief Gravity, MRF, and wall-heat mathematical source contributions.

#include "SF_sourceContribution.h"

#include "core/system/SF_systemContribution.h"

#include <algorithm>
#include <array>
#include <stdexcept>

namespace SF::Physics::SourceContribution {
namespace {

struct SourceTarget {
    FDM::SourceKind kind;
    const char* symbol;
    bool energy;
};

constexpr std::array<SourceTarget,3> targets{{
    {FDM::SourceKind::Gravity,"gravity",false},
    {FDM::SourceKind::MRF,"MRF",false},
    {FDM::SourceKind::WallHeat,"wallHeat",true}
}};

} // namespace

void contribute(
        System::SystemContribution& system,
        const FDM::SourceConfig& config) {
    for (FDM::SourceKind kind : config.enabled) {
        const auto item = std::find_if(
            targets.begin(),targets.end(),
            [kind](const SourceTarget& value) { return value.kind == kind; });
        if (item == targets.end()) {
            throw std::runtime_error(
                "No equation contribution is registered for SourceKind.");
        }
        system.recordContribution(
            std::string("model.")+item->symbol,
            std::string(item->symbol)+" equation source contribution");
        auto source=System::FormulaExpr::op("source",
            {System::FormulaExpr::symbol(item->symbol)},item->symbol);
        system.extendMathematics(item->energy
            ? std::vector<std::string>{"energy","enthalpy","E_ENTHALPY"}
            : std::vector<std::string>{"momentum"},std::move(source));
    }
}

} // namespace SF::Physics::SourceContribution
