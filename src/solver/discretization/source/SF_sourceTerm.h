/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

/// @file SF_sourceTerm.h
/// @brief 源项离散算子 (explicit source terms).
///
/// 源项写入 Field::Source，后续由 Math::spatialResidual 以 −Source 计入残差。
/// 调用方式:
/// @code
///   SourceTerm::Sp(field, sourceConfig);
/// @endcode

#include "SF_field.h"
#include "core/residual/SF_residual.h"
#include "SF_config.h"
#include "methods/numerics/structured/SF_structured.h"
#include "SF_valueTypes.h"
#include "gravity/SF_gravity.h"
#include "heat/SF_wallHeatSource.h"
#include "mrf/SF_mrf.h"

#include <algorithm>
#include <array>
#include <stdexcept>
#include <string>

namespace SF {

/// @brief 源项算子命名空间 (显式源项)。
///        使用 "SourceTerm" 与 physics 下各具体 `Source::*` 模型区分。
namespace SourceTerm {

/// @brief Check whether a source setting applies to a cell.
/// @param field Field containing zone/set membership.
/// @param zone Zone name; empty or `"all"` applies globally.
/// @param i Cell i index.
/// @param j Cell j index.
/// @param k Cell k index.
/// @return True when the cell should receive the source contribution.
inline bool appliesToZone(const Field& field, const std::string& zone, int i, int j, int k) {
    if (zone.empty() || zone == "all" || zone == "All" || zone == "ALL") return true;
    return Source::MRF::appliesToZone(field, zone, i, j, k);
}

/// @brief Add gravity source terms from explicit config.
/// @param field Field whose source array is updated in-place.
/// @param settings Zone-scoped gravity accelerations.
inline void addGravity(Field& field, Residual& residual,
                       const std::vector<ZoneVectorSetting>& settings) {
    if (settings.empty()) {
        throw std::runtime_error(
            "Gravity source is enabled but no acceleration setting was loaded.");
    }

    Math::forFluidInterior(field, [&](int i, int j, int k) {
        for (const auto& setting : settings) {
            if (appliesToZone(field, setting.zone, i, j, k)) {
                Source::Gravity::addTranslation(field, residual, i, j, k, setting.value);
            }
        }
    });
}

/// @brief Add MRF/rotating-frame source terms from explicit config.
/// @param field Field whose source array is updated in-place.
/// @param settings Zone-scoped rotating-frame settings.
inline void addMRF(Field& field, Residual& residual,
                   const std::vector<RotatingSetting>& settings) {
    if (settings.empty()) {
        throw std::runtime_error(
            "MRF source is enabled but no rotating setting was loaded.");
    }

    Math::forFluidInterior(field, [&](int i, int j, int k) {
        for (const auto& setting : settings) {
            if (!appliesToZone(field, setting.zone, i, j, k)) continue;

            Vector3 frameVelocity = setting.hasVelocity ? setting.velocity : Vector3();
            Source::Rotating::addRotatingFrame(field, residual, i, j, k,
                                               setting.center,
                                               Source::MRF::angularVelocity(setting),
                                               frameVelocity);
        }
    });
}

/// @brief Add wall heat sources from explicit config.
/// @param field Field whose energy source is updated in-place.
/// @param settings Patch-scoped wall heat flux settings.
inline void addWallHeat(Field& field, Residual& residual,
                        const std::vector<WallHeatSetting>& settings) {
    Source::WallHeat::addSource(field, residual, settings);
}

using ContributionAssembler = void(*)(
    Field&, Residual&, const FDM::SourceConfig&);

struct RegisteredContribution {
    FDM::SourceKind kind;
    const char* name;
    ContributionAssembler assemble;
};

inline const std::array<RegisteredContribution,3>& contributions() {
    static const std::array<RegisteredContribution,3> registry{{
        {FDM::SourceKind::Gravity,"gravity",
         [](Field& field, Residual& residual,
            const FDM::SourceConfig& config) {
             addGravity(field,residual,config.gravity);
         }},
        {FDM::SourceKind::MRF,"MRF",
         [](Field& field, Residual& residual,
            const FDM::SourceConfig& config) {
             addMRF(field,residual,config.rotating);
         }},
        {FDM::SourceKind::WallHeat,"wallHeat",
         [](Field& field, Residual& residual,
            const FDM::SourceConfig& config) {
             addWallHeat(field,residual,config.wallHeat);
         }}
    }};
    return registry;
}

/// @brief Assemble source terms from explicit solver config.
///
/// This overload is preferred by the Equation layer because it receives all source data
/// as a value object and does not read source-selection globals.
///
/// @param field Field whose source array is cleared then updated in-place.
/// @param config Explicit source configuration.
inline void Sp(Field& field, Residual& residual,
               const FDM::SourceConfig& parameters,
               const std::vector<FDM::SourceKind>& boundSources) {
    residual.clearSource();

    for (FDM::SourceKind kind : boundSources) {
        const auto found = std::find_if(
            contributions().begin(),contributions().end(),
            [kind](const RegisteredContribution& item) {
                return item.kind == kind;
            });
        if (found == contributions().end()) {
            throw std::runtime_error(
                "No density equation contribution is registered for SourceKind.");
        }
        found->assemble(field,residual,parameters);
    }
}

} // namespace SourceTerm

} // namespace SF
