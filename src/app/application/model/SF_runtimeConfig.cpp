/// @file SF_runtimeConfig.cpp
/// @brief 构造 mesh 与 IBM 所需的不可变 runtime configuration。
///
/// Data flow:
///   parsed case/model configuration
///       -> ghost-layer and boundary-set requirements
///       -> mesh/IBM initialization configuration
///
/// 本文件不建立 MPI communicator、不执行数值算法，也不写运行报告。

#include "app/application/model/SF_runtimeConfig.h"

#include "SF_multiphase.h"

#include <algorithm>
#include <vector>

namespace SF::Application::Runtime {
namespace {

void appendSetName(std::vector<std::string>& names, const std::string& name) {
    if (name.empty()) return;
    if (std::find(names.begin(), names.end(), name) == names.end()) {
        names.push_back(name);
    }
}

template <typename T>
void appendSetNames(
        std::vector<std::string>& names,
        const std::vector<BCSetting<T>>& settings) {
    for (const auto& setting : settings) appendSetName(names, setting.name);
}

void appendMultiPhaseSetNames(
        std::vector<std::string>& names,
        const Physics::Multiphase::MultiPhaseConfig& config) {
    for (const auto& assignment : config.setPhases) {
        appendSetName(names, assignment.setName);
    }
    for (const auto& plane : config.phiPlaneInitializers) {
        appendSetName(names, plane.setName);
    }
    for (const auto& sphere : config.phiSphereInitializers) {
        appendSetName(names, sphere.setName);
    }
    for (const auto& condition : config.phiInitialConditions) {
        appendSetName(names, condition.setName);
    }
    for (const auto& condition : config.phiBoundaryConditions) {
        appendSetName(names, condition.setName);
    }
    appendSetNames(names, config.alpha.initialConditions);
    appendSetNames(names, config.alpha.boundaryConditions);
    appendSetNames(names, config.temperature.initialConditions);
    appendSetNames(names, config.temperature.boundaryConditions);
}

} // namespace

MeshRuntimeConfig makeMeshRuntimeConfig(
        const CaseConfig& caseConfig,
        const FDM::SolverConfig& solverConfig,
        int requiredTermHaloWidth) {
    MeshRuntimeConfig config;
    config.initialConditions = solverConfig.initial;
    config.idealGasGamma = solverConfig.numerics.idealGasGamma;
    config.idealGasConstant = solverConfig.numerics.idealGasConstant;
    config.convectionScheme = solverConfig.numerics.recipes.convection
        ? FDM::toString(solverConfig.numerics.recipes.convection->id())
        : "unbound";
    config.ilwOrder = solverConfig.numerics.ilwOrder;
    config.requiredGhostLayers = std::max(
        requiredTermHaloWidth,
        FDM::requiredGhostLayersForILW(config.ilwOrder));
    config.parallelEnabled = caseConfig.parallel.enabled;
    config.parallelProcessCount = caseConfig.parallel.processCount;
    config.automaticPartition = caseConfig.parallel.automaticPartition;
    config.partitionSplit = caseConfig.parallel.partitionSplit;
    config.haloTolerance = caseConfig.parallel.haloTolerance;

    appendSetNames(config.configuredSetNames, solverConfig.initial.velocity);
    appendSetNames(config.configuredSetNames, solverConfig.initial.pressure);
    appendSetNames(config.configuredSetNames, solverConfig.initial.density);
    appendSetNames(config.configuredSetNames, solverConfig.initial.temperature);
    appendSetNames(config.configuredSetNames,
                   solverConfig.turbulence.scalars.kInitial);
    appendSetNames(config.configuredSetNames,
                   solverConfig.turbulence.scalars.epsilonInitial);
    appendSetNames(config.configuredSetNames,
                   solverConfig.turbulence.scalars.omegaInitial);
    appendSetNames(config.configuredSetNames, solverConfig.boundaries.velocity);
    appendSetNames(config.configuredSetNames,
                   solverConfig.boundaries.energyFromPressure);
    appendSetNames(config.configuredSetNames, solverConfig.boundaries.density);
    for (const auto& setting : solverConfig.boundaries.thermal) {
        appendSetName(config.configuredSetNames, setting.name);
    }
    appendSetNames(config.configuredSetNames,
                   solverConfig.turbulence.scalars.kBoundary);
    appendSetNames(config.configuredSetNames,
                   solverConfig.turbulence.scalars.epsilonBoundary);
    appendSetNames(config.configuredSetNames,
                   solverConfig.turbulence.scalars.omegaBoundary);
    appendSetNames(config.physicalBoundaryNames,
                   solverConfig.boundaries.density);
    appendSetNames(config.physicalBoundaryNames,
                   solverConfig.boundaries.velocity);
    appendSetNames(config.physicalBoundaryNames,
                   solverConfig.boundaries.energyFromPressure);
    for (const auto& setting : solverConfig.boundaries.thermal) {
        appendSetName(config.physicalBoundaryNames, setting.name);
    }
    appendSetNames(config.physicalBoundaryNames,
                   solverConfig.turbulence.scalars.kBoundary);
    appendSetNames(config.physicalBoundaryNames,
                   solverConfig.turbulence.scalars.epsilonBoundary);
    appendSetNames(config.physicalBoundaryNames,
                   solverConfig.turbulence.scalars.omegaBoundary);
    if (caseConfig.multiPhaseEnabled) {
        appendMultiPhaseSetNames(
            config.configuredSetNames, caseConfig.multiPhase);
    }
    return config;
}

IBM::IBMRuntimeConfig makeIBMRuntimeConfig(
        const FDM::SolverConfig& solverConfig,
        int requiredTermHaloWidth) {
    IBM::IBMRuntimeConfig config;
    config.enabled = solverConfig.ibm.enabled;
    config.requiredGhostLayers = std::max(
        requiredTermHaloWidth,
        FDM::requiredGhostLayersForILW(solverConfig.numerics.ilwOrder));
    config.ilwEnabled =
        solverConfig.numerics.ibmBoundary == FDM::IBMBoundaryScheme::ILW;
    config.ilwAccuracyOrder = solverConfig.numerics.ilwOrder;
    config.gamma = solverConfig.numerics.idealGasGamma;
    config.normalAngleDegrees = solverConfig.ibm.normalAngleDegrees;
    config.normalSearchMinLayers = solverConfig.ibm.normalSearchMinLayers;
    config.normalSearchMaxLayers = solverConfig.ibm.normalSearchMaxLayers;
    config.normalSearchTargetCandidates =
        solverConfig.ibm.normalSearchTargetCandidates;
    config.normalSearchKeepSamples = solverConfig.ibm.normalSearchKeepSamples;
    config.method = solverConfig.ibm.method;
    config.forcing = solverConfig.ibm.forcing;
    return config;
}

} // namespace SF::Application::Runtime
