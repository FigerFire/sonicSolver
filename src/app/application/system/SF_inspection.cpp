/// @file SF_inspection.cpp
/// @brief 构造稳定的 case inspection snapshot。

#include "SF_inspection.h"

#include "SF_algorithmDescriptor.h"
#include "SF_systemBuilder.h"
#include "SF_systemValidator.h"
#include "app/application/model/SF_model.h"
#include "app/application/model/SF_runtimeConfig.h"

#include <stdexcept>

namespace SF::Application {
namespace {

/// @brief 从 case 的 physics/equation template 直接解析物理状态族。
///
/// 数据依赖方向是 CaseConfig -> ResolvedSimulationSystem；执行计划不参与
/// physics/equation template 的解析。
System::PhysicsTemplateKind resolvePhysics(const CaseConfig& config) {
    if (!config.multiPhaseEnabled || !config.multiPhase.enabled) {
        return System::PhysicsTemplateKind::SingleFluid;
    }

    const std::string& type = config.multiPhase.type;
    if (Physics::Multiphase::isHomogeneousType(type)) {
        return System::PhysicsTemplateKind::HomogeneousMixture;
    }
    if (Physics::Multiphase::isEulerianEulerianType(type)) {
        return System::PhysicsTemplateKind::EulerianEulerian;
    }
    if (Physics::Multiphase::isLevelSetType(type)) {
        return System::PhysicsTemplateKind::OneFluidInterface;
    }

    // legacy mixture-like 及未显式匹配的多相模型，现有 builder 统一落到
    // mixture-like 输入仍组合 homogeneous-mixture 物理状态族。
    return System::PhysicsTemplateKind::HomogeneousMixture;
}

std::string turbulenceName(FDM::TurbulenceModelKind model) {
    switch (model) {
        case FDM::TurbulenceModelKind::None: return "none";
        case FDM::TurbulenceModelKind::kEpsilon: return "kEpsilon";
        case FDM::TurbulenceModelKind::kOmegaSST: return "kOmegaSST";
        case FDM::TurbulenceModelKind::Smagorinsky: return "Smagorinsky";
        case FDM::TurbulenceModelKind::DNS: return "DNS";
    }
    throw std::runtime_error("Unknown turbulence model.");
}

} // namespace

CaseInspection inspectCase(const CaseConfig& config) {
    CaseInspection result;
    result.config = config;
    result.mesh = Runtime::makeMeshRuntimeConfig(config, config.solver);
    result.ibm = Runtime::makeIBMRuntimeConfig(config.solver);

    // 数据依赖方向：CaseConfig -> ResolvedSimulationSystem。
    // 物理状态族直接从 case 的 physics/template 输入解析，不再由
    // legacy WorkflowPlan.kind 决定。
    FDM::ImmersedAlgorithmDescriptor immersedDescriptor;
    System::BuildRequest systemRequest;
    systemRequest.templateOrigin = resolvePhysics(config);
    systemRequest.phaseNames =
        config.multiPhase.eulerianEulerian.phaseNames;
    systemRequest.levelSet =
        config.multiPhaseEnabled && config.multiPhase.enabled
        && Physics::Multiphase::isLevelSetType(config.multiPhase.type);
    systemRequest.homogeneousThermodynamics =
        config.multiPhaseEnabled && config.multiPhase.enabled
        && Physics::Multiphase::isHomogeneousType(config.multiPhase.type);
    systemRequest.legacyMixture =
        config.multiPhaseEnabled && config.multiPhase.enabled
        && Physics::Multiphase::isMixtureType(config.multiPhase.type);
    systemRequest.phaseChange = config.multiPhase.phaseChange.enabled;
    systemRequest.transportedLegacyAlpha =
        systemRequest.legacyMixture
        && config.multiPhase.alpha.transportEnabled;
    systemRequest.interfaceGhostFluid = systemRequest.levelSet
        && Physics::Multiphase::normalizeModelType(
            config.multiPhase.levelSet.surfaceTensionModel) == "ghostfluid";
    systemRequest.turbulence = config.solver.turbulence.enabled;
    systemRequest.turbulenceModel = turbulenceName(config.solver.turbulence.model);
    systemRequest.turbulencePhaseNames =
        config.solver.turbulence.phaseNames;
    systemRequest.parallel = config.parallel.enabled;
    // Runtime 已提供 ConstraintGlobalDof 的 canonical owner、连续 row layout
    // 和 owner-to-consumer COPY。它可用于 Peskin/FTS/fractional DLM 的分布式
    // surface/body constraint；pressure monolithic KKT 是否可进入仍由 algorithm
    // stage 的更细粒度 fail-fast 负责，不能在 case inspection 中一概拒绝。
    systemRequest.distributedLinearSystemAvailable = true;
    systemRequest.constraintGlobalDofAvailable = true;
    systemRequest.composition = config.composition;
    if (result.ibm.enabled) {
        immersedDescriptor = result.ibm.method == FDM::IBMMethod::Ghost
            ? IBM::Descriptor::ghostCell()
            : IBM::Descriptor::variational(result.ibm.forcing);
        systemRequest.immersed = &immersedDescriptor;
    }
    result.system = System::build(config.solver, systemRequest);
    System::validate(result.system);

    return result;
}

CaseInspection inspectCase(const std::string& path) {
    return inspectCase(ModelLoader::read(path));
}

} // namespace SF::Application
