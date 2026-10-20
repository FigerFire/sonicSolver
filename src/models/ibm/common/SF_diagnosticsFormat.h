#pragma once

/// @file SF_diagnosticsFormat.h
/// @brief Ghost/ILW setup 诊断文本格式化。

#include "SF_ibmConfig.h"
#include "method/ghost/SF_weightBuilder.h"

#include <string>

namespace SF::IBM::Common {

std::string formatILWConfiguration(const IBMRuntimeConfig& config);
std::string formatILWPlans(const GhostIBM::WeightBuilder& weights);
std::string formatPreprocessingTiming(
    const GhostIBM::WeightBuilder& weights, double stlSeconds);

} // namespace SF::IBM::Common
