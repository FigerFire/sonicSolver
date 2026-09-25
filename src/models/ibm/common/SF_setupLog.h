#pragma once

/// @file SF_setupLog.h
/// @brief IBM setup 的结构化终端日志输出。

#include "SF_ibmConfig.h"
#include "SF_immersedSystem.h"
#include "geoProcessing/SF_STLGeometry.h"
#include "method/ghost/SF_weightBuilder.h"

namespace SF::IBM::Common {

void logGeometryDiagnostics(const GeoProcessing::STLGeometry& geometry);
void logGhostSetup(const IBMRuntimeConfig& config,
                   const GhostIBM::WeightBuilder& weights,
                   double stlSeconds);
void logForcingSetup(const FDM::ImmersedMethodSelection& selection);

} // namespace SF::IBM::Common
