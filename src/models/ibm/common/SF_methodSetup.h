#pragma once

/// @file SF_methodSetup.h
/// @brief 根据单一 method 选择完成 ghost 或 forcing 方法装配。

#include "SF_field.h"
#include "SF_ibmConfig.h"
#include "SF_immersedSystem.h"
#include "geoProcessing/SF_STLGeometry.h"
#include "method/ghost/SF_weightBuilder.h"
#include "method/SF_method.h"

namespace SF::IBM::Common {

FDM::ImmersedAlgorithmDescriptor setupMethod(
    Field& field,
    const GeoProcessing::STLGeometry& geometry,
    const IBMRuntimeConfig& config,
    double stlSeconds,
    GhostIBM::WeightBuilder& ghostWeights,
    Forcing::ImmersedForcingSystem& forcing);

} // namespace SF::IBM::Common
