/// @file SF_diagnosticsFormat.cpp
/// @brief Ghost/ILW setup 诊断文本格式化实现。

#include "common/SF_diagnosticsFormat.h"

#include "common/SF_timingFormat.h"

#include <sstream>

namespace SF::IBM::Common {

std::string formatILWConfiguration(const IBMRuntimeConfig& config) {
    std::ostringstream out;
    out << "accuracy/Taylor=" << config.ilwAccuracyOrder << "/"
        << config.requestedTaylorOrder()
        << ", normalAngle=" << config.normalAngleDegrees
        << " deg, alignmentThreshold="
        << config.normalAlignmentThreshold()
        << ", searchLayers=" << config.normalSearchMinLayers
        << ".." << config.normalSearchMaxLayers
        << ", targetCandidates=" << config.normalSearchTargetCandidates
        << ", keepSamples=" << config.normalSearchKeepSamples;
    return out.str();
}

std::string formatILWPlans(const GhostIBM::WeightBuilder& weights) {
    const auto& stats = weights.ilwStats();
    std::ostringstream out;
    out << "order>=1 " << stats.built << "/" << stats.candidates
        << ", ghostCandidates=" << stats.ghostCandidates
        << ", solidCandidates=" << stats.solidCandidates
        << ", reused=" << stats.reusedPlans
        << ", directBuilt=" << stats.directBuiltPlans
        << ", order>=2=" << stats.built - stats.maxOrderLessThanTwo
        << ", order>=4=" << stats.built - stats.maxOrderLessThanFour
        << ", order>=6=" << stats.built - stats.maxOrderLessThanSix
        << ", order>=8=" << stats.built - stats.maxOrderLessThanEight
        << ", fluidSampleFail=" << stats.insufficientFluidSamples
        << ", curvatureFail=" << stats.curvatureFailures
        << ", fitFail=" << stats.fitFailures;
    return out.str();
}

std::string formatPreprocessingTiming(
        const GhostIBM::WeightBuilder& weights, double stlSeconds) {
    const auto& timing = weights.timing();
    const auto& stats = weights.ilwStats();
    std::ostringstream out;
    out << "stlLoadBVH=" << formatSeconds(stlSeconds)
        << ", sdfClassify=" << formatSeconds(timing.sdfSeconds)
        << ", ghostLayer=" << formatSeconds(timing.ghostLayerSeconds)
        << ", flagGeometry=" << formatSeconds(timing.classifySeconds)
        << ", ilwSampleCache=" << formatSeconds(timing.ilwSampleCacheSeconds)
        << " (cells=" << timing.ilwSampleCells
        << ", fluid=" << formatSeconds(timing.ilwFluidSampleSeconds)
        << ", normal=" << formatSeconds(timing.ilwNormalSampleSeconds) << ")"
        << ", ilwPlan=" << formatSeconds(timing.ilwPlanSeconds)
        << " (geometry=" << formatSeconds(stats.planGeometrySeconds)
        << ", fluidSamples=" << formatSeconds(stats.planSampleSeconds)
        << ", spacing=" << formatSeconds(stats.planSpacingSeconds)
        << ", probe=" << formatSeconds(stats.planProbeSeconds)
        << ", wallFit=" << formatSeconds(stats.planWallFitSeconds)
        << ", curvature=" << formatSeconds(stats.planCurvatureSeconds)
        << ", higherFit=" << formatSeconds(stats.planHigherFitSeconds) << ")"
        << ", donorWeights=" << formatSeconds(timing.donorWeightSeconds)
        << ", buildTotal=" << formatSeconds(timing.totalSeconds)
        << ", setupTotal=" << formatSeconds(stlSeconds + timing.totalSeconds);
    return out.str();
}

} // namespace SF::IBM::Common
