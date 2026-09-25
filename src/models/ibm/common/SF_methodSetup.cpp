/// @file SF_methodSetup.cpp
/// @brief IBM ghost/forcing 方法装配实现。

#include "common/SF_methodSetup.h"

#include "common/SF_methodSelection.h"
#include "common/SF_setupLog.h"
#include "descriptor/SF_algorithmDescriptor.h"
#include "core/interfaces/SF_log.h"

namespace SF::IBM::Common {

FDM::ImmersedAlgorithmDescriptor setupMethod(
        Field& field,
        const GeoProcessing::STLGeometry& geometry,
        const IBMRuntimeConfig& config,
        double stlSeconds,
        GhostIBM::WeightBuilder& ghostWeights,
        Forcing::ImmersedForcingSystem& forcing) {
    logGeometryDiagnostics(geometry);
    if (config.method == FDM::IBMMethod::VariationalForcing) {
        field.clearCellFlags(FLUID_CELL);
        forcing.configure(geometry, config);
        logForcingSetup(selectMethod(config));
        return forcing.algorithmDescriptor();
    }

    broadcast("IBM preprocessing: ", "start");
    ghostWeights.build(field, geometry, config);
    logGhostSetup(config, ghostWeights, stlSeconds);
    return Descriptor::ghostCell();
}

} // namespace SF::IBM::Common
