/// @file SF_setupLog.cpp
/// @brief IBM setup 的结构化终端日志输出实现。

#include "common/SF_setupLog.h"

#include "common/SF_diagnosticsFormat.h"
#include "core/interfaces/SF_log.h"

namespace SF::IBM::Common {

void logGeometryDiagnostics(const GeoProcessing::STLGeometry& geometry) {
    if (geometry.watertight()) return;
    const auto& report = geometry.watertightReport();
    broadcast("IBM warning: ",
              "STL is not watertight. SDF signs may be unreliable.");
    broadcast("IBM boundary edges: ", report.boundaryEdges);
    broadcast("IBM non-manifold edges: ", report.nonManifoldEdges);
}

void logGhostSetup(const IBMRuntimeConfig& config,
                   const GhostIBM::WeightBuilder& weights,
                   double stlSeconds) {
    const auto& counts = weights.counts();
    broadcast("IBM fluid cells: ", counts.fluid);
    broadcast("IBM ghost cells: ", counts.ghost);
    broadcast("IBM solid cells: ", counts.solid);
    const auto& stats = weights.ilwStats();
    if (stats.candidates > 0 || config.ilwEnabled) {
        if (config.ilwEnabled) {
            broadcast("IBM ILW configuration: ",
                      formatILWConfiguration(config));
        }
        broadcast("IBM ILW plans: ", formatILWPlans(weights));
    }
    broadcast("IBM preprocessing timing: ",
              formatPreprocessingTiming(weights, stlSeconds));
    broadcast("IBM preprocessing: ", "complete");
    broadcast("IBM setup: ", "ghost-cell method ready");
}

void logForcingSetup(const FDM::ImmersedMethodSelection& selection) {
    broadcast(
        "IBM setup: ",
        selection.enforcement == FDM::IBMEnforcement::MonolithicKKT
        ? "surface variational KKT ready; extended fluid domain remains "
          "canonical FLUID_CELL state"
        : "variational forcing strategy ready; extended fluid domain "
          "remains canonical FLUID_CELL state");
}

} // namespace SF::IBM::Common
