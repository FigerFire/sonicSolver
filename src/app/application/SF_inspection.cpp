#include "solver/system/SF_systemPrinter.h"
#include <sstream>
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
#include "models/physics/gravity/SF_accelerationProvider.h"
#include "models/physics/mrf/SF_frameProvider.h"
#include "models/physics/heat/SF_wallFluxProvider.h"
#include "models/physics/multiphase/SF_multiphaseSystemContribution.h"
#include "models/ibm/SF_ibmSystemContribution.h"

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

    if (Physics::Multiphase::isThermalType(type)) return System::PhysicsTemplateKind::SingleFluid;
    if (Physics::Multiphase::isMixtureType(type)) return System::PhysicsTemplateKind::HomogeneousMixture;
    throw std::runtime_error("Unknown physical system contribution: "+type);
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
    systemRequest.configuredSources=true;
    systemRequest.userContributions=config.userContributions;
    // Explicit WHAT decides whether the default fluid bundle is present;
    // adding an unrelated user source must not change storage selection.
    if(config.composition.declared || !config.userContributions.empty()) {
        systemRequest.includeDefaultFluidPreset=std::any_of(config.composition.equations.begin(),config.composition.equations.end(),[](const auto& name) {
            return name=="Continuity" || name=="Momentum" || name=="Energy" || name=="PressureConstraint";
        });
    }
    systemRequest.templateOrigin = resolvePhysics(config);
    systemRequest.authoredExecution=config.authoredExecution;
    systemRequest.authoredNumerics=config.authoredNumerics;
    systemRequest.phaseNames =
        config.multiPhase.eulerianEulerian.phaseNames;
    systemRequest.referencePhase = config.multiPhase.eulerianEulerian.referencePhase;
    const bool levelSet = config.multiPhaseEnabled && config.multiPhase.enabled
        && Physics::Multiphase::isLevelSetType(config.multiPhase.type);
    if (config.multiPhaseEnabled && config.multiPhase.enabled && Physics::Multiphase::isThermalType(config.multiPhase.type)
        && (config.parallel.enabled || config.meshFiles.size()!=1))
        throw std::runtime_error("Unsupported thermal transport: serial single patch boundary/halo binding required.");
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
    systemRequest.homogeneousTurbulenceUnsupported = config.solver.turbulence.enabled
        && config.solver.turbulence.family != FDM::TurbulenceFamily::DNS;
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
            config.solver.pressure,true);
    }
    systemRequest.parallel = config.parallel.enabled;
    // host backend 能力来自真实链接的 backend，而不是 SystemBuilder 的假设。
    systemRequest.capabilities = detectBuildCapabilities();
    systemRequest.composition = config.composition;
    systemRequest.modelContributions=ModelLoader::modelDeclarations(config);
    if (result.ibm.enabled) {
        immersedDescriptor = result.ibm.method == FDM::IBMMethod::Ghost
            ? IBM::Descriptor::ghostCell() : IBM::Descriptor::variational(result.ibm.forcing);
        immersedDescriptor.wallClosure=config.solver.ibm.wallClosure;
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
    // Direct named provider storage declares its native initialization/BC backing.
    // Reconstructed phase primaries and specialized ports keep their owning binder contracts.
    if(config.modelDescription)for(const auto& symbol:result.system.executableSystem.state.symbols()) {
        if(symbol.storageBinding!=System::StorageBinding::ProviderDistributed
            || symbol.storageKey!=symbol.id || symbol.role!=System::StateRole::Transported)continue;
        const auto& fields=config.modelDescription->fields;
        const auto field=std::find_if(fields.begin(),fields.end(),[&](const auto& value){return value.symbol==symbol.id;});
        if(field==fields.end())throw std::runtime_error("Required native STATE lacks registered field: "+symbol.id);
        if(symbol.initializationRequired && field->initial.is_null())
            throw std::runtime_error("Required native STATE initialization missing: "+symbol.id+" ("+symbol.origin.source+")");
        if(symbol.boundaryRequired && field->boundaries.empty())
            throw std::runtime_error("Required native STATE boundary missing: "+symbol.id+" ("+symbol.origin.source+")");
    }
    result.mesh = Runtime::makeMeshRuntimeConfig(
        config,config.solver,result.system.numericalSystem.requiredHaloWidth);
    const auto* pressureInputs=System::requiredPressureInputs(System::programOf(result.system));
    result.mesh.initializeConservativeState=System::requiredPackedComponents(result.system.executableSystem.state,
        result.system.solvePlan.compiledProgram,"conservative")>0
        || System::requiredPhaseInputs(System::programOf(result.system))!=nullptr;
    if(config.modelDescription)for(const auto& field:config.modelDescription->fields) {
        if(!result.system.executableSystem.state.contains(field.symbol))continue;
        if((result.mesh.initializeConservativeState || pressureInputs)
            && result.system.executableSystem.state.at(field.symbol).storageBinding!=System::StorageBinding::NamedDistributed)continue;
        // A declared scalar default initializes the entire geometry, unlike a fluid set-only IC.
        if(field.initial.is_number())result.mesh.configuredSetNames.push_back("all");
        for(const auto& item:field.sets.items())result.mesh.configuredSetNames.push_back(item.key());
        for(const auto& item:field.boundaries.items()) {
            result.mesh.configuredSetNames.push_back(item.key());result.mesh.physicalBoundaryNames.push_back(item.key());
        }
    }
    for(const auto& symbol:result.system.executableSystem.state.symbols())if(symbol.fieldInput) {
        for(const auto& value:symbol.fieldInput->initialConditions)result.mesh.configuredSetNames.push_back(value.name);
        for(const auto& value:symbol.fieldInput->boundaries) {
            result.mesh.configuredSetNames.push_back(value.name);result.mesh.physicalBoundaryNames.push_back(value.name);
        }
    }
    result.ibm = Runtime::makeIBMRuntimeConfig(
        config.solver,result.system.numericalSystem.requiredHaloWidth);

    return result;
}

CaseInspection inspectCase(const std::string& path) {
    return inspectCase(ModelLoader::read(path));
}

std::string describeCase(const CaseInspection& inspection) {
    std::ostringstream out;out << "BUILT-IN COMPOSITION / SELECTION ORIGINS\n"
        << "  optional model: absent=disabled; present=enabled unless enabled:false; laminar adds no turbulence transport.\n";
    if(inspection.config.modelDescription)for(const auto& [target,source]:inspection.config.modelDescription->selectionOrigins)
        out << "  " << target << " <- " << source << "\n";
    const auto& config=inspection.config.solver;
    out << "  actual turbulence: " << (config.turbulence.enabled ? "active" : "inactive (laminar/no turbulence transport)")
        << "; model=" << FDM::toString(config.turbulence.model) << "\n";
    out << "  actual source contributors: gravity=" << config.sources.gravity.size()
        << " MRF=" << config.sources.rotating.size() << " wallHeat=" << config.sources.wallHeat.size() << "\n";
    out << "  inactive model parameters/field files may remain authored; compiled consumers determine resources.\n";
    if(inspection.config.modelDescription) {
        out<<"REGISTERED FIELD BOUNDARY INPUTS (runtime resources validated on run)\n";
        for(const auto& field:inspection.config.modelDescription->fields)if(inspection.system.executableSystem.state.contains(field.symbol))
            for(const auto& item:field.boundaries.items())out<<"  "<<field.symbol<<" / "<<item.key()<<" : "<<item.value().dump()<<"\n";
    }
    return out.str()+System::describe(inspection.system)+IBM::renderExplain(inspection.ibmExplain);
}

} // namespace SF::Application
