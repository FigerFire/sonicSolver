#include "core/system/SF_operationIds.h"
/// @file SF_systemBuilder.cpp
/// @brief Compile built-in equations and neutral model/user contribution values.

#include "SF_systemBuilder.h"
#include "SF_couplingStatus.h"
#include "SF_equationContribution.h"
#include "SF_solvePlan.h"
#include "SF_transformation.h"
#include "SF_numericalCompiler.h"
#include "SF_presets.h"
#include "SF_singleFluidPreset.h"
#include "SF_pressureCoupling.h"
#include "SF_builtinState.h"
#include "SF_executionComposition.h"
#include "SF_eulerianRelations.h"
#include "SF_providerResolver.h"

#include "SF_config.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace SF::System {
namespace {


bool hasRawEquationPrefix(
        const RawEquationSystem& system, std::string_view prefix) {
    return std::any_of(system.legacyEquations.begin(),system.legacyEquations.end(),
        [&](const EquationDescriptor& value) {
            return value.id.compare(0,prefix.size(),prefix) == 0;
        });
}

bool usesEquation(
        const EquationCompositionConfig& composition,
        std::string_view name) {
    return std::find(composition.equations.begin(),composition.equations.end(),name)
        != composition.equations.end();
}

bool hasOperation(const RuntimeReport& runtime, std::string_view id) {
    return std::find(
        runtime.requiredOperations.begin(),runtime.requiredOperations.end(),id)
        != runtime.requiredOperations.end();
}

bool hasExecutableEquation(
        const ResolvedSimulationSystem& system, std::string_view id) {
    return std::any_of(
        system.executableSystem.legacyEquations.begin(),
        system.executableSystem.legacyEquations.end(),
        [&](const EquationDescriptor& equation) { return equation.id == id; });
}

void addState(
        SystemCompositionBuilder& system,
        std::string id,
        std::string name,
        int components = 1,
        StateRole role = StateRole::Primary,
        StorageBinding binding = StorageBinding::SpecializedExecutor,
        std::string storageKey = {},
        int componentOffset = 0,
        std::string nameSpace = {}) {
    StateSymbol unknown;
    unknown.id = std::move(id);
    unknown.name = std::move(name);
    unknown.components = components;
    unknown.shape = components == 1 ? ValueShape::Scalar : ValueShape::Vector;
    unknown.role = role;
    unknown.storageBinding = binding;
    unknown.storageKey = std::move(storageKey);
    unknown.componentOffset = componentOffset;
    unknown.nameSpace = std::move(nameSpace);
    system.addState(std::move(unknown));
}

void addEquation(
        SystemCompositionBuilder& system,
        EquationDescriptor descriptor,
        SF::Equation::Definition definition) {
    system.addEquation(std::move(descriptor),std::move(definition));
}

void requireProvider(
        ResolvedSimulationSystem& system, std::string id,
        std::string responsibility) {
    const auto found = std::find_if(
        system.runtime.providerRequirements.begin(),system.runtime.providerRequirements.end(),
        [&](const ProviderRequirement& value) { return value.id == id; });
    if (found == system.runtime.providerRequirements.end()) {
        system.runtime.providerRequirements.push_back(
            {std::move(id),std::move(responsibility)});
    }
}

void requireRuntimeService(
        ResolvedSimulationSystem& system, std::string id,
        std::string reason) {
    const auto found = std::find_if(
        system.runtime.runtimeServiceRequirements.begin(),
        system.runtime.runtimeServiceRequirements.end(),
        [&](const RuntimeServiceRequirement& value) { return value.id == id; });
    if (found == system.runtime.runtimeServiceRequirements.end()) {
        system.runtime.runtimeServiceRequirements.push_back(
            {std::move(id),std::move(reason)});
    }
}

void deriveExecutionComposition(
        ResolvedSimulationSystem& system, const BuildRequest& request) {
    const bool eulerianOperations = std::any_of(
        system.runtime.operationBindings.begin(),
        system.runtime.operationBindings.end(),
        [](const ResolvedOperationBinding& binding) {
            return binding.provider == "flow.eulerian-pressure";
        });
    for (const auto& binding : system.runtime.operationBindings) {
        if (binding.status == BindingStatus::Resolved) {
            requireProvider(system,binding.provider,
                            "execute compiled operation '"+binding.operation+"'");
        }
    }

    if (hasExecutableEquation(system,"E_LEGACY_ALPHA")) {
        requireProvider(system,"equation.legacy-mixture",
                        "bind transported mixture state and RHS contribution");
    }
    if (hasExecutableEquation(system,"E_PHAScontinuity")) {
        requireProvider(system,"thermodynamics.homogeneous",
                        "bind homogeneous primary state and equation set");
    } else if (requiresProvider(system,"flow.conservative")) {
        requireProvider(system,"thermodynamics.single-fluid",
                        "bind the conservative single-fluid equation set");
    }
    if (std::any_of(
            system.runtime.requirements.begin(),system.runtime.requirements.end(),
            [](const ExecutionRequirement& requirement) {
                return requirement.name == "PhaseChangeExecution"
                    && requirement.required;
            })) {
        requireProvider(system,"equation.phase-change",
                        "bind phase-change equation contribution");
    }
    for (const auto* contributions:
         {&request.modelContributions,&request.userContributions})
        for (const auto& contribution:*contributions)
            for (const auto& required:contribution.providerRequirements)
                requireProvider(system,required.id,required.responsibility);

    requireRuntimeService(system,"execution.contract",
                          "honor halo, canonical COPY, and SUM contracts");
    if (request.parallel) {
        requireRuntimeService(system,"mpi.distributed",
                              "execute distributed ownership contracts");
    }
    const bool needsHypre = eulerianOperations
        || std::any_of(
            system.executableSystem.operations.begin(),
            system.executableSystem.operations.end(),
            [&](const ExecutableOperation& operation) {
                return hasOperation(system.runtime.report,operation.operation)
                    && std::find(operation.requirements.begin(),
                                 operation.requirements.end(),
                                 OperationCapability::PressureLinearSolve)
                        != operation.requirements.end();
            })
        || hasOperation(system.runtime.report,OpIds::IbmKktSolve);
    if (needsHypre) {
        requireProvider(system,"linear.hypre",
                        "bind pressure or monolithic distributed linear solve");
        requireRuntimeService(system,"mpi.initialized",
                              "HYPRE requires an initialized MPI context");
    }
}

void validateSemanticComposition(const EquationCompositionConfig& composition) {
    if (!composition.declared) return;
    if (composition.equations.empty()) {
        throw std::runtime_error(
            "equations.yaml is required when semantic composition is declared.");
    }
    const bool continuity = usesEquation(composition,"Continuity");
    const bool momentum = usesEquation(composition,"Momentum");
    const bool energy = usesEquation(composition,"Energy");
    if (!continuity || !momentum) {
        throw std::runtime_error(
            "Current single-fluid composition requires Continuity and Momentum.");
    }
    const auto providers = builtinCompositionProviders();
    if (composition.thermoDynamics.equationOfState.empty()) {
        throw std::runtime_error(
            "Active fluid equations require models/thermoDynamics.yaml with equationOfState.");
    }
    if (!providers.hasEquationOfState(composition.thermoDynamics.equationOfState)) {
        throw std::runtime_error("No registered equation-of-state provider for '"
            +composition.thermoDynamics.equationOfState+"'.");
    }
    if (energy && composition.thermoDynamics.thermo.empty()) {
        throw std::runtime_error(
            "Energy equation requires caloric thermodynamics, but no thermo model is registered.");
    }
    if (energy && !providers.hasCaloricThermo(composition.thermoDynamics.thermo)) {
        throw std::runtime_error("No registered caloric-thermo provider for '"
            +composition.thermoDynamics.thermo+"'.");
    }
    if (!composition.thermoDynamics.transport.empty()
        && !providers.hasTransport(composition.thermoDynamics.transport))
        throw std::runtime_error("No registered transport provider for '"+composition.thermoDynamics.transport+"'.");
    if (!composition.stateDeclared || composition.solutionVariables.empty())
        throw std::runtime_error("Native composition requires explicit solution STATE (stateRegistry.use).");
    if (composition.algorithm.empty() || !providers.hasAlgorithm(composition.algorithm))
        throw std::runtime_error("Native composition requires explicit HOW (algorithmRegistry).");

}

} // namespace

ResolvedSimulationSystem build(
        const FDM::SolverConfig& config,
        const BuildRequest& request) {
    if (!request.unsupportedCompositionReason.empty())
        throw std::runtime_error(request.unsupportedCompositionReason);
    ResolvedSimulationSystem result;
    ExecutionProgram executionProgram;
    std::vector<NumericalBinding> equationMethods;
    result.classification.templateOrigin = request.templateOrigin;
    result.timeRecipe = config.numerics.timeRecipe;
    validateSemanticComposition(request.composition);
    if (request.composition.declared) {
        const auto& algorithm=request.composition.algorithm;
        if ((algorithm=="Explicit" && request.coupling)
            || (algorithm!="Explicit" && (!request.coupling
                || algorithm!=FDM::toString(request.coupling->presetKind))))
            throw std::runtime_error("Explicit HOW selection and contributed coupling contract disagree.");
    }
    const bool legacyFluidDiffusion = config.numerics.viscousEnabled
        || request.modelRequiresDiffusion || request.legacyMixture;

    std::vector<TransformationDescriptor> transformationRequests;
    SystemCompositionBuilder builtin(
        result.rawSystem,transformationRequests,result.executionPolicies,
        {request.templateOrigin==PhysicsTemplateKind::EulerianEulerian ? OriginKind::Model : OriginKind::BuiltinPreset,
         request.templateOrigin==PhysicsTemplateKind::EulerianEulerian ? "Eulerian PhaseSystem" : "default system"});

    const auto& selection=request.composition;
    result.rawSystem.state.setSelectionOrigin(selection.stateSelectionOrigin);
    if (selection.declared && request.templateOrigin!=PhysicsTemplateKind::SingleFluid)
        throw std::runtime_error("Native single-fluid equation/STATE composition has no current multiphase provider contract.");
    if (request.templateOrigin!=PhysicsTemplateKind::EulerianEulerian) {
        if (!selection.stateDeclared || selection.solutionVariables.empty())
            throw std::runtime_error("Single-fluid equations require explicit solution STATE selection.");
        const auto includes=[&](const char* id) {
            return std::find(selection.solutionVariables.begin(),selection.solutionVariables.end(),id)!=selection.solutionVariables.end();
        };
        if (includes("U") && includes("rhoU"))
            throw std::runtime_error("Ambiguous Momentum target: solution STATE selects both U and rhoU; provide an explicit equation-target contract.");
        if (!includes("U") && !includes("rhoU"))
            throw std::runtime_error("Momentum has no selected solution STATE target (U or rhoU).");
        const BuiltinStateCatalog catalog;
        for (const auto& id:selection.solutionVariables)
            if (catalog.contains(id)) {
                if (result.rawSystem.state.isSolution(id))
                    throw std::runtime_error("Duplicate solution STATE: "+id);
                builtin.addState(catalog.solution(id));
                result.rawSystem.state.selectSolution(id);
            }
        if (selection.declared && selection.thermoDynamics.equationOfState=="rhoConst"
            && !includes("rho")) {
            auto density=catalog.at("rho");
            density.name="constant density closure";
            density.storageBinding=StorageBinding::SpecializedExecutor;
            density.storageKey="rhoConst";
            density.constantValue=selection.thermoDynamics.constantDensity;
            density.restartEligible=false;
            builtin.addState(std::move(density));
        }
        if (const auto reason=validateSolutionProviderContract(result.rawSystem.state,selection); !reason.empty())
            throw std::runtime_error(reason);
        if (selection.declared && includes("U")) {
            if (!includes("p")) throw std::runtime_error("Pressure constraint requires selected solution p; HOW cannot add a solution unknown.");
            Compose::addConstantDensityFluid(builtin,selection,legacyFluidDiffusion);
        } else if (request.pressureConstraint || selection.compatibilityPressureConstraint
            || usesEquation(selection,"PressureConstraint")) {
            Compose::addPressureConstraintFluid(builtin,PressureConstraintSpec{legacyFluidDiffusion},selection);
        } else {
            Preset::installSingleFluid(builtin,SingleFluidPresetSpec{
                request.singleFluidPreset ? request.singleFluidPreset->diffusion : legacyFluidDiffusion},selection);
        }
    } else {
        Compose::addEulerianEulerianTemplate(builtin,result,request.phaseNames,
            request.referencePhase.empty() && !request.phaseNames.empty() ? request.phaseNames.front() : request.referencePhase);
    }

    if (request.composition.declared) {
        const auto& eos = request.composition.thermoDynamics.equationOfState;
        if (eos == "rhoConst") {
            result.classification.densityBehavior = "constant";
            result.classification.thermodynamicCompressibility = "zero";
            result.rawSystem.closures.push_back("rho = rho0 (rhoConst)");
        } else if (eos == "perfectGas") {
            result.classification.densityBehavior = "variable";
            result.classification.thermodynamicCompressibility = "available";
            result.rawSystem.closures.push_back("p = p(rho,T) (perfectGas)");
        }
    }

    if (result.classification.densityBehavior.empty()) {
        result.classification.densityBehavior = "legacy-configured";
        result.classification.thermodynamicCompressibility = "legacy-configured";
    }

    SystemCompositionBuilder models(
        result.rawSystem,transformationRequests,result.executionPolicies,
        {OriginKind::Model,"configured models"});
    if (request.homogeneousThermodynamics) {
        models.recordContribution(
            "model.homogeneousMultiphase","homogeneous multiphase equations");
        addState(models,"phaseMassAux","homogeneous phase mass auxiliary",1,
                   StateRole::Transported,
                   StorageBinding::SpecializedExecutor,{},0,"homogeneous");
        addEquation(models,{
            "E_PHAScontinuity","homogeneous phase-mass transport",
            "conservation",{"phaseMassAux"}},
            SF::Equation::named("E_PHAScontinuity",
                SF::Equation::ddt({"phaseMassAux"})
                    + SF::Equation::div({"phaseMassFlux"})
                    == SF::Equation::Symbol{"phaseMassSources"}));
    } else if (request.legacyMixture) {
        models.recordContribution(
            "model.legacyMixture","legacy mixture equations");
        addState(models,"alphaAux","legacy transported volume fraction",1,
                   StateRole::Transported,
                   StorageBinding::SpecializedExecutor,{},0,"legacyMultiphase");
        addEquation(models,{
            "E_LEGACY_ALPHA","legacy volume-fraction transport",
            "conservation",{"alphaAux"}},
            SF::Equation::named("E_LEGACY_ALPHA",
                SF::Equation::ddt({"alphaAux"})
                    + SF::Equation::div({"alphaFlux"})
                    == SF::Equation::Symbol{"zero"}));
    }
    for (const SystemContribution& contribution:request.modelContributions)
        models.applyContribution(contribution);
    // §15 precedence：用户修改必须排在 builtin defaults 与 model
    // contributions 之后、formulation/transformation 之前。当前没有 typed
    // lowering 的动作会显式失败，绝不静默覆盖。
    SystemCompositionBuilder users(
        result.rawSystem,transformationRequests,result.executionPolicies,
        {OriginKind::User,"user C++ contribution"});
    for (const SystemContribution& contribution:request.userContributions)
        users.applyContribution(contribution);
    for (const SystemModification& modification : request.userModifications) {
        users.applyModification(modification);
    }

    const BuiltinStateCatalog stateCatalog;
    for (const auto* module:{&builtin,&models,&users})
        for (const auto& symbol:module->requiredStates)
            stateCatalog.require(result.rawSystem.state,symbol);

    if (request.templateOrigin!=PhysicsTemplateKind::EulerianEulerian)
    for (const auto& id:selection.solutionVariables)
        if (!result.rawSystem.state.isSolution(id)) result.rawSystem.state.selectSolution(id);

    TransformerRegistry transformers;
    const auto momentum=std::find_if(builtin.execution.begin(),builtin.execution.end(),[](const auto& node) {
        return node.kind==ExecutionKind::EquationCall && node.step.equation=="momentum";
    });
    transformers.registerTransformer(makePressureConstraintTransformer(
        momentum==builtin.execution.end()?std::string{}:momentum->step.target.symbol));
    transformers.registerTransformer(makeSharedPressureTransformer());
    transformers.registerTransformer(makeImmersedConstraintTransformer());
    // REGISTERED MODEL：压力耦合 preset。它的生效条件来自 resolved raw
    // equation/constraint structure，而不是 density/pressure 标签；状态必须
    // 显式报告，禁止静默忽略。
    if (request.coupling) {
        result.coupling = contributePressureCoupling(
            models,*request.coupling,result.rawSystem);
    }
    // phase-1 是否真的贡献了 policy；resolution 之后 status 可能变成
    // Unsupported（缺 formulation operation），但 plan 仍必须真实描述该
    // schedule，provider resolver 才能报告缺失。
    const bool couplingPolicyContributed = isCouplingActive(result.coupling);
    const bool monolithicPolicy = std::any_of(
        result.executionPolicies.begin(),result.executionPolicies.end(),
        [](const LegacyExecutionPolicy& policy) {
            return policy.kind==LegacyExecutionPolicyKind::MonolithicKKT;
        });
    if (monolithicPolicy
        && couplingPolicyContributed) {
        // DLM/KKT 是完全隐式约束耦合：它取代 segregated pressure policy，
        // 不是与之并存。preset 仍然必须被显式报告为未生效。
        result.executionPolicies.erase(
            std::remove_if(
                result.executionPolicies.begin(),result.executionPolicies.end(),
                [](const LegacyExecutionPolicy& policy) {
                    return policy.id == kPressureScheduleId;
                }),
            result.executionPolicies.end());
        result.coupling.status = CouplingStatus::Inactive;
        result.coupling.reason =
            "superseded by the monolithic KKT/DLM constraint contribution";
    }
    result.executableSystem = TransformationPipeline::apply(
        result.rawSystem,transformationRequests,transformers,
        result.executionPolicies,result.transformations);
    if (hasConstraint(result.executableSystem,"C_SHARED_PRESSURE"))
        projectEulerianBackendDefinitions(result.executableSystem);
    for (auto* module:{&builtin,&models,&users}) {
        executionProgram.root.children.insert(executionProgram.root.children.end(),
            module->execution.begin(),module->execution.end());
        equationMethods.insert(equationMethods.end(),module->numerics.begin(),module->numerics.end());
        executionProgram.legacyEntries.insert(executionProgram.legacyEntries.end(),
            module->legacyExecution.begin(),module->legacyExecution.end());
    }
    const auto legacyFragments=composeContributions(result,request,executionProgram,equationMethods);
    requireTargetStates(result.executableSystem.state,executionProgram.root);
    result.numericalSelection=selectNumerics(config,request,result.executableSystem,std::move(equationMethods));
    // Source composition is complete; the compiler receives immutable WHAT/HOW/WHICH.
    NumericalCompiler::compileSystem(result,executionProgram,result.numericalSelection,legacyFragments,
        request.additionalExecutionContributions);
    for (const auto& provider:result.numericalSystem.providerRequirements)
        requireProvider(result,provider,"execute a compiled mathematical term binding");
    if (!result.solvePlan.compiledProgram.steps.empty() && result.coupling.status==CouplingStatus::Active)
        result.coupling.derivedOperations=
            SolvePlanner::requiredOperations(result.solvePlan);
    // A declared equation without HOW lowering is an explicit capability gap;
    // it must never inherit the Flow timestep merely for report ownership.
    for (const auto& equation:result.executableSystem.legacyEquations) {
        const bool owned=std::any_of(result.solvePlan.blocks.begin(),
            result.solvePlan.blocks.end(),[&](const SolveBlock& block) {
                return std::find(block.equations.begin(),block.equations.end(),
                                 equation.id)!=block.equations.end();
            });
        if (!owned) {
            result.runtime.report.status=RuntimeStatus::Unsupported;
            result.runtime.report.reason += " Equation '"+equation.id
                +"' has no compiled HOW execution step.";
        }
    }

    const bool nativeTurbulence=std::any_of(result.runtime.operationBindings.begin(),
        result.runtime.operationBindings.end(),[](const auto& binding) {
            return binding.operation==OpIds::TurbulenceAdvance;
        });
    result.runtime.requirements.push_back({"SingleFluidTurbulenceScope",nativeTurbulence,
        !nativeTurbulence || (!request.parallel && request.templateOrigin==PhysicsTemplateKind::SingleFluid
            && std::none_of(result.rawSystem.contributions.begin(),result.rawSystem.contributions.end(),
                [](const auto& contribution) { return contribution.id.rfind("model.ibm.",0)==0; })),
        "native RAS transport is single-fluid, serial, single patch; pressure/IBM/MPI providers are unavailable"});
    result.runtime.requirements.push_back({
        "EulerianGlobalDof",true,true,"canonical cell ownership"});
    // 共享面的 canonical flux identity 需要真实的 rank-to-rank 通信；
    // serial stub 不构成并行能力。
    result.runtime.requirements.push_back({
        "CanonicalFace",request.parallel,
        !request.parallel || request.capabilities.mpi,
        "one owner computes each shared face flux"});
    // Backend availability 是 host 事实（BuildCapabilities），case 只回答
    // “是否需要”。编译器把两者相比，得到 required/available/reason 三者一致
    // 的 RuntimeRequirement，并在 validation 阶段 fail-fast。
    const bool constraintDof = request.parallel && std::any_of(
        result.rawSystem.state.symbols().begin(),result.rawSystem.state.symbols().end(),
        [](const StateSymbol& unknown) {
            return unknown.location==VariableLocation::BodyConstraint
                || unknown.location==VariableLocation::SurfaceConstraint;
        });
    result.runtime.requirements.push_back({
        "ConstraintGlobalDof",constraintDof,
        !constraintDof || request.capabilities.canonicalConstraintDof,
        "unique owner for Lambda/lambda and sparse J/S edges"});
    const bool distributedSolve = monolithicPolicy;
    result.runtime.requirements.push_back({
        "DistributedLinearSystem",distributedSolve,
        !distributedSolve || request.capabilities.distributedLinearSystem,
        "GlobalDofId is mapped to backend rows outside equation assembly"});
    const bool homogeneousUnsupported = request.homogeneousThermodynamics
        && (!result.realization.conservativeTransportedMass
            || request.homogeneousTurbulenceUnsupported);
    result.runtime.requirements.push_back({
        "HomogeneousEquationExecution",request.homogeneousThermodynamics,
        !homogeneousUnsupported,
        "current homogeneous FluidStateModel execution requires density formulation "
        "without transported turbulence"});
    const bool eulerianExecution =
        request.templateOrigin == PhysicsTemplateKind::EulerianEulerian;
    const bool eulerianDtPolicy =
        std::isfinite(config.numerics.maxDeltaT)
        && config.numerics.maxDeltaT > 0.0;
    result.runtime.requirements.push_back({
        "EulerianEquationExecution",eulerianExecution,
        !eulerianExecution || eulerianDtPolicy,
        "current Eulerian equation executor requires a shared-pressure "
        "constraint formulation and a finite positive maxDeltaT"});
    result.runtime.requirements.push_back({
        "PhaseChangeExecution",request.phaseChange,
        request.phaseChangeExecutionAvailable,
        "phase change is implemented for homogeneous or Eulerian equation systems"});
    result.runtime.requirements.push_back({
        "CanonicalScalarInterfaceFlux",
        request.transportedLegacyAlpha && request.parallel,
        !(request.transportedLegacyAlpha && request.parallel),
        "legacy transported alpha on coupled patches is not implemented"});
    deriveExecutionComposition(result,request);
    return result;
}

bool hasUnknown(const ResolvedSimulationSystem& system, std::string_view id) {
    return hasUnknown(system.executableSystem,id);
}

bool hasUnknown(const RawEquationSystem& system, std::string_view id) {
    return std::any_of(system.state.symbols().begin(),system.state.symbols().end(),
        [&](const StateSymbol& value) { return value.id == id; });
}

bool hasUnknown(const ExecutableEquationSystem& system, std::string_view id) {
    return std::any_of(system.state.symbols().begin(),system.state.symbols().end(),
        [&](const StateSymbol& value) { return value.id == id; });
}

bool hasEquation(const ResolvedSimulationSystem& system, std::string_view id) {
    return hasEquation(system.executableSystem,id);
}

bool hasEquation(const RawEquationSystem& system, std::string_view id) {
    return system.registry.contains(id);
}

bool hasEquation(const ExecutableEquationSystem& system, std::string_view id) {
    return system.registry.contains(id);
}

const SF::Equation::Definition& equationDefinition(
        const ResolvedSimulationSystem& system, std::string_view id) {
    return system.executableSystem.legacyDefinitions.at(std::string(id));
}

bool hasEquationPrefix(
        const ExecutableEquationSystem& system, std::string_view prefix) {
    return std::any_of(
        system.legacyEquations.begin(),
        system.legacyEquations.end(),
        [&](const EquationDescriptor& value) {
            return value.id.compare(0, prefix.size(), prefix) == 0;
        });
}

bool hasEquationPrefix(
        const ResolvedSimulationSystem& system, std::string_view prefix) {
    return hasEquationPrefix(system.executableSystem,prefix);
}

bool hasConstraint(const ResolvedSimulationSystem& system, std::string_view id) {
    return hasConstraint(system.executableSystem,id);
}

bool hasConstraint(const RawEquationSystem& system, std::string_view id) {
    return std::any_of(system.constraints.begin(),system.constraints.end(),
        [&](const ConstraintDescriptor& value) { return value.id == id; });
}

bool hasConstraint(
        const ExecutableEquationSystem& system, std::string_view id) {
    return std::any_of(system.constraints.begin(),system.constraints.end(),
        [&](const ConstraintDescriptor& value) { return value.id == id; });
}

bool requiresCapability(
        const RuntimeRequirements& requirements, std::string_view name) {
    return std::any_of(
        requirements.requirements.begin(), requirements.requirements.end(),
        [&](const ExecutionRequirement& value) {
            return value.name == name && value.required;
        });
}

bool requiresProvider(
        const RuntimeRequirements& requirements, std::string_view id) {
    return std::any_of(
        requirements.providerRequirements.begin(),
        requirements.providerRequirements.end(),
        [&](const ProviderRequirement& value) { return value.id == id; });
}

bool requiresRuntimeService(
        const RuntimeRequirements& requirements, std::string_view id) {
    return std::any_of(
        requirements.runtimeServiceRequirements.begin(),
        requirements.runtimeServiceRequirements.end(),
        [&](const RuntimeServiceRequirement& value) { return value.id == id; });
}

bool requiresCapability(
        const ResolvedSimulationSystem& system, std::string_view name) {
    return requiresCapability(system.runtime, name);
}

bool requiresProvider(
        const ResolvedSimulationSystem& system, std::string_view id) {
    return requiresProvider(system.runtime, id);
}

bool requiresRuntimeService(
        const ResolvedSimulationSystem& system, std::string_view id) {
    return requiresRuntimeService(system.runtime, id);
}

} // namespace SF::System
