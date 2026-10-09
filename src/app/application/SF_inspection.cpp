/// @file SF_inspection.cpp
/// @brief 构造稳定的 case inspection snapshot。

#include "SF_inspection.h"

#include "SF_algorithmDescriptor.h"
#include "SF_hostCapabilities.h"
#include "SF_pressureCoupling.h"
#include "SF_systemBuilder.h"
#include "SF_systemValidator.h"
#include "SF_termProviderCatalog.h"
#include "app/application/model/SF_model.h"
#include "app/application/model/SF_runtimeConfig.h"
#include "models/physics/SF_sourceContribution.h"
#include "models/physics/gravity/SF_accelerationProvider.h"
#include "models/physics/mrf/SF_frameProvider.h"
#include "models/physics/heat/SF_wallFluxProvider.h"
#include "models/physics/interfaceModel/levelSet/SF_levelSetSystemContribution.h"
#include "models/physics/multiphase/SF_multiphaseSystemContribution.h"
#include "models/ibm/SF_ibmSystemContribution.h"
#include "models/turbulence/SF_turbulenceSystemContribution.h"

#include <algorithm>
#include <stdexcept>

namespace SF::Application {
namespace {

/// 把原始 CaseConfig 转换成一份已经解析、校验、补全了运行需求的 CaseInspection 快照
// case 文件
//    ↓
// CaseConfig
//    ↓
// inspectCase()
//    ├─ 判断物理体系
//    ├─ 解析 IBM
//    ├─ 构造 BuildRequest
//    ├─ build ResolvedSimulationSystem
//    ├─ validate
//    ├─ 根据数值系统反推出 halo width
//    └─ 生成最终 Mesh / IBM runtime config
//    ↓
// CaseInspection
/// 数据依赖方向是 CaseConfig -> ResolvedSimulationSystem；执行计划不参与
/// physics/equation template 的解析。

// resolvePhysics根据 case 的 multiphase 配置，判断当前物理状态属于哪个 physics template
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

// @brief 将 turbulence model 枚举转换为字符串，便于在 CaseInspection 中使用。
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
    // IBM algorithm selection is independent of halo depth. The final runtime
    // config is rebuilt below from the compiled numerical requirement.
    result.ibm = Runtime::makeIBMRuntimeConfig(config.solver,0);
    if (result.ibm.enabled && result.ibm.wallClosure==FDM::ImmersedWallClosure::StationaryNoSlipAdiabatic) {
        if (config.parallel.enabled || config.meshFiles.size()>1
            || resolvePhysics(config)!=System::PhysicsTemplateKind::SingleFluid)
            throw std::runtime_error("stationaryNoSlipAdiabatic currently requires a serial single-fluid patch; MPI/phase wall ports are Unsupported.");
        if (result.ibm.method!=FDM::IBMMethod::Ghost || result.ibm.ilwEnabled)
            throw std::runtime_error("stationaryNoSlipAdiabatic is the explicit quadratic Ghost closure; Euler ILW compatibility is not its numerical method. Set ILW off explicitly.");
        if (config.solver.numerics.viscous!=FDM::ViscousScheme::Central2)
            throw std::runtime_error("stationaryNoSlipAdiabatic currently validates only Central2 viscous gradients.");
        if (!config.solver.numerics.viscousEnabled)
            throw std::runtime_error("stationaryNoSlipAdiabatic requires active molecular viscous diffusion.");
        if (config.solver.turbulence.enabled && config.solver.turbulence.model!=FDM::TurbulenceModelKind::kOmegaSST
            && config.solver.turbulence.family!=FDM::TurbulenceFamily::DNS)
            throw std::runtime_error("stationaryNoSlipAdiabatic turbulence port currently implements only SST.");
    }

    // 数据依赖方向：CaseConfig -> ResolvedSimulationSystem。
    // 物理状态族直接从 case 的 physics/template 输入解析，不再由
    // legacy WorkflowPlan.kind 决定。
    FDM::ImmersedAlgorithmDescriptor immersedDescriptor;
    System::BuildRequest systemRequest;
    systemRequest.templateOrigin = resolvePhysics(config);
    systemRequest.authoredExecution=config.authoredExecution;
    systemRequest.authoredNumerics=config.authoredNumerics;
    systemRequest.phaseNames =
        config.multiPhase.eulerianEulerian.phaseNames;
    systemRequest.referencePhase = config.multiPhase.eulerianEulerian.referencePhase;
    const bool levelSet = config.multiPhaseEnabled && config.multiPhase.enabled
        && Physics::Multiphase::isLevelSetType(config.multiPhase.type);
    if (levelSet) {
        System::SystemContribution contribution;
        Physics::InterfaceModels::LevelSetContribution::contribute(
            contribution,{Physics::Multiphase::normalizeModelType(
                config.multiPhase.levelSet.surfaceTensionModel)=="ghostfluid",
                config.multiPhase.levelSet.advectionOrder,
                config.multiPhase.levelSet.reinitializationOrder,
                config.multiPhase.levelSet.reinitializationSteps,
                config.multiPhase.levelSet.pseudoTimeStep,
                config.multiPhase.levelSet.wenoEpsilon,
                config.multiPhase.levelSet.wenoPower,
                config.multiPhase.levelSet.signSmoothingFactor,
                Physics::Multiphase::normalizeModelType(config.multiPhase.levelSet.surfaceTensionModel)=="csf",
                config.multiPhase.levelSet.surfaceTension,config.multiPhase.levelSet.interfaceThickness});
        systemRequest.modelContributions.push_back(std::move(contribution));
    }
    systemRequest.homogeneousThermodynamics =
        config.multiPhaseEnabled && config.multiPhase.enabled
        && Physics::Multiphase::isHomogeneousType(config.multiPhase.type);
    systemRequest.legacyMixture =
        config.multiPhaseEnabled && config.multiPhase.enabled
        && Physics::Multiphase::isMixtureType(config.multiPhase.type);
    if (systemRequest.homogeneousThermodynamics || systemRequest.legacyMixture) {
        System::SystemContribution contribution;
        if (systemRequest.homogeneousThermodynamics)
            Physics::Multiphase::contributeHomogeneous(contribution,config.multiPhase);
        else Physics::Multiphase::contributeMixture(contribution,config.multiPhase);
        systemRequest.modelContributions.push_back(std::move(contribution));
    }
    systemRequest.phaseChange = config.multiPhase.phaseChange.enabled;
    systemRequest.transportedLegacyAlpha =
        systemRequest.legacyMixture
        && config.multiPhase.alpha.transportEnabled;
    if (config.solver.turbulence.enabled) {
        systemRequest.homogeneousTurbulenceUnsupported =
            config.solver.turbulence.family != FDM::TurbulenceFamily::DNS;
        System::SystemContribution contribution;
        Turbulence::contribute(contribution,Turbulence::SystemContributionSpec{
            turbulenceName(config.solver.turbulence.model),
            config.solver.turbulence.phaseNames,
            systemRequest.templateOrigin
                == System::PhysicsTemplateKind::EulerianEulerian,
            systemRequest.phaseNames});
        systemRequest.modelContributions.push_back(std::move(contribution));
    }
    systemRequest.modelRequiresDiffusion =
        config.solver.turbulence.enabled || systemRequest.legacyMixture
        || levelSet;
    const bool fluidDiffusion = config.solver.numerics.viscousEnabled
        || systemRequest.modelRequiresDiffusion;
    if (systemRequest.templateOrigin==System::PhysicsTemplateKind::SingleFluid
        && !config.composition.declared) {
        if (config.composition.compatibilityPressureConstraint)
            systemRequest.pressureConstraint=System::PressureConstraintSpec{fluidDiffusion};
        else
            systemRequest.singleFluidPreset=System::SingleFluidPresetSpec{fluidDiffusion};
    }
    if (config.pressureCouplingDeclared) {
        // 显式注册的 coupling preset：是否生效由 resolved equation system
        // 决定，状态在 explain/validation 中显式报告。
        systemRequest.coupling = System::couplingRequestFrom(
            config.solver.pressure.coupling,true);
    }
    systemRequest.parallel = config.parallel.enabled;
    // host backend 能力来自真实链接的 backend，而不是 SystemBuilder 的假设。
    systemRequest.capabilities = detectBuildCapabilities();
    systemRequest.composition = config.composition;
    {
        System::SystemContribution contribution;
        Physics::SourceContribution::contribute(
            contribution,config.solver.sources);
        systemRequest.modelContributions.insert(
            systemRequest.modelContributions.begin(),std::move(contribution));
    }
    if (result.ibm.enabled) {
        immersedDescriptor = result.ibm.method == FDM::IBMMethod::Ghost
            ? IBM::Descriptor::ghostCell()
            : IBM::Descriptor::variational(result.ibm.forcing);
        immersedDescriptor.wallClosure=config.solver.ibm.wallClosure;
        System::SystemContribution contribution;
        IBM::SystemContribution::contribute(contribution,immersedDescriptor);
        systemRequest.modelContributions.push_back(std::move(contribution));
    }
    systemRequest.additionalExecutionContributions = result.ibm.enabled
        || config.solver.turbulence.enabled || levelSet
        || systemRequest.legacyMixture
        || systemRequest.homogeneousThermodynamics
        || systemRequest.phaseChange
        || systemRequest.transportedLegacyAlpha;
    systemRequest.phaseChangeExecutionAvailable =
        !(systemRequest.phaseChange
          && (levelSet || systemRequest.legacyMixture));
    for (auto descriptor:Physics::Gravity::termProviders(
             config.solver.sources.gravity))
        systemRequest.termProviders.push_back(std::move(descriptor));
    for (auto descriptor:Physics::MRF::termProviders(
             config.solver.sources.rotating))
        systemRequest.termProviders.push_back(std::move(descriptor));
    systemRequest.termProviders.push_back(
        Physics::WallHeat::termProvider(config.solver.sources.wallHeat));
    result.ibmExplain = IBM::explainSnapshot(immersedDescriptor);
    result.system = System::build(config.solver, systemRequest);
    System::validate(result.system);
    result.mesh = Runtime::makeMeshRuntimeConfig(
        config,config.solver,result.system.numericalSystem.requiredHaloWidth);
    result.ibm = Runtime::makeIBMRuntimeConfig(
        config.solver,result.system.numericalSystem.requiredHaloWidth);

    return result;
}

CaseInspection inspectCase(const std::string& path) {
    return inspectCase(ModelLoader::read(path));
}

} // namespace SF::Application
