/// @file SF_systemBuilder.cpp
/// @brief 显式声明单流体、多相、湍流和 IBM 对数学系统的贡献。

#include "SF_systemBuilder.h"
#include "SF_couplingStatus.h"
#include "SF_equationContribution.h"
#include "SF_couplingStatus.h"
#include "SF_solvePlan.h"
#include "SF_couplingStatus.h"
#include "SF_transformation.h"
#include "SF_couplingStatus.h"
#include "SF_numericalCompiler.h"
#include "SF_providerResolver.h"
#include "SF_couplingStatus.h"
#include "SF_presets.h"
#include "SF_singleFluidPreset.h"
#include "SF_couplingStatus.h"
#include "SF_pressureCoupling.h"
#include "SF_couplingStatus.h"
#include "SF_stateRealization.h"
#include "SF_couplingStatus.h"
#include "models/physics/SF_sourceContribution.h"
#include "models/ibm/SF_ibmSystemContribution.h"
#include "SF_immersedSystem.h"

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
    return std::any_of(system.equations.begin(),system.equations.end(),
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
        system.executableSystem.equations.begin(),
        system.executableSystem.equations.end(),
        [&](const EquationDescriptor& equation) { return equation.id == id; });
}

bool hasExecutableEquationPrefix(
        const ExecutableEquationSystem& system, std::string_view prefix) {
    return std::any_of(
        system.equations.begin(),
        system.equations.end(),
        [&](const EquationDescriptor& equation) {
            return equation.id.compare(0,prefix.size(),prefix) == 0;
        });
}

void addUnknown(
        SystemCompositionBuilder& system,
        std::string id,
        std::string name,
        int components = 1,
        UnknownRole role = UnknownRole::Primary,
        StorageBinding binding = StorageBinding::SpecializedExecutor,
        std::string storageKey = {},
        int componentOffset = 0,
        std::string nameSpace = {}) {
    UnknownDescriptor unknown;
    unknown.id = std::move(id);
    unknown.name = std::move(name);
    unknown.components = components;
    unknown.shape = components == 1 ? ValueShape::Scalar : ValueShape::Vector;
    unknown.role = role;
    unknown.storageBinding = binding;
    unknown.storageKey = std::move(storageKey);
    unknown.componentOffset = componentOffset;
    unknown.nameSpace = std::move(nameSpace);
    system.addUnknown(std::move(unknown));
}

void addEquation(
        SystemCompositionBuilder& system,
        EquationDescriptor descriptor,
        Equation::Definition definition) {
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
    if (hasExecutableEquation(system,"E_LEVEL_SET")) {
        requireProvider(system,"equation.level-set",
                        "bind interface state, transport, jump, and curvature");
    }
    if (hasExecutableEquation(system,"E_PHASE_MASS")) {
        requireProvider(system,"thermodynamics.homogeneous",
                        "bind homogeneous primary state and equation set");
    } else if (!eulerianOperations) {
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
    const bool turbulenceClosure = std::any_of(
        system.executableSystem.closures.begin(),
        system.executableSystem.closures.end(),
        [](const std::string& closure) {
            return closure.compare(0,4,"mu_t") == 0
                || closure.find("turbulence") != std::string::npos;
        });
    if (turbulenceClosure) {
        requireProvider(system,"closure.turbulence",
                        "bind turbulence closure and transported equations");
        if (hasExecutableEquationPrefix(system.executableSystem,"E_TURB_")) {
            requireProvider(system,"equation.turbulence-transport",
                            "bind transported turbulence state and equations");
        }
    }
    if (request.immersed) {
        if (request.immersed->enforcement == FDM::IBMEnforcement::GhostCell) {
            requireProvider(system,"ibm.boundary",
                            "bind ghost/ILW boundary closure");
        }
        if (hasOperation(system.runtime.report,"ibm.constraint.project")
            || hasOperation(system.runtime.report,"ibm.kkt.solve")) {
            requireProvider(system,"ibm.constraint",
                            "bind immersed constraint state and operations");
        }
    }

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
        || hasOperation(system.runtime.report,"ibm.kkt.solve");
    if (needsHypre) {
        requireProvider(system,"linear.hypre",
                        "bind pressure or monolithic distributed linear solve");
        requireRuntimeService(system,"mpi.initialized",
                              "HYPRE requires an initialized MPI context");
    }
}

FDM::PressureCouplingPreset pressurePreset(const std::string& kind) {
    if (kind == "SIMPLE") return FDM::PressureCouplingPreset::SIMPLE;
    if (kind == "PISO") return FDM::PressureCouplingPreset::PISO;
    if (kind == "PIMPLE") return FDM::PressureCouplingPreset::PIMPLE;
    throw std::runtime_error(
        "A pressure-constraint composition requires SIMPLE, PISO, or PIMPLE.");
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
    if (composition.thermoDynamics.equationOfState == "rhoConst") {
        if (energy) throw std::runtime_error(
            "rhoConst composition currently supports Momentum + Continuity; "
            "select an explicit Energy/Enthalpy equation before adding caloric closure.");
        if (!providers.hasTransport(composition.thermoDynamics.transport)) {
            throw std::runtime_error(
                "rhoConst Momentum requires transport; no transport model is registered.");
        }
        // rhoConst 方程族不选择 coupling preset；preset 由 BuildRequest 显式
        // 注册并由 resolved constraint structure 匹配。这里只要求输入没有
        // 声明一个无法解析的 preset 名字。
        if (!composition.algorithm.empty()) {
            (void)pressurePreset(composition.algorithm);
        }
    } else if (composition.thermoDynamics.equationOfState == "perfectGas") {
        if (!energy) throw std::runtime_error(
            "perfectGas single-fluid execution currently requires Energy; "
            "a non-caloric variable-density equation pack is unsupported.");
        if (!composition.algorithm.empty() && composition.algorithm != "Explicit") {
            throw std::runtime_error(
                "The selected pressure algorithm requires an incompressibility constraint; "
                "the perfectGas conservative system has none.");
        }
    }
}

} // namespace

ResolvedSimulationSystem build(
        const FDM::SolverConfig& config,
        const BuildRequest& request) {
    ResolvedSimulationSystem result;
    result.classification.templateOrigin = request.templateOrigin;
    result.timeRecipe = config.numerics.timeRecipe;
    validateSemanticComposition(request.composition);
    const bool legacyFluidDiffusion = config.numerics.viscousEnabled
        || request.turbulence || request.legacyMixture || request.levelSet;

    std::vector<TransformationDescriptor> transformationRequests;
    SystemCompositionBuilder builtin(
        result.rawSystem,transformationRequests,result.executionPolicies,
        {OriginKind::BuiltinPreset,"default system"});

    if (request.singleFluidPreset) {
        result.classification.templateOrigin = PhysicsTemplateKind::SingleFluid;
        Preset::installSingleFluid(
            builtin,*request.singleFluidPreset);
    } else if (request.composition.declared) {
        if (request.templateOrigin != PhysicsTemplateKind::SingleFluid) {
            throw std::runtime_error(
                "Semantic single-fluid equation composition cannot be combined "
                "with the current multiphase template.");
        }
        // Equations choose the raw mathematical system. EOS contributes a
        // closure and may simplify that system; it never selects a solver
        // family or swaps in a different equation pack.
        if (usesEquation(request.composition,"Energy")) {
            throw std::runtime_error(
                "Density single-fluid composition requires the explicit "
                "single-fluid preset selection.");
        }
        Compose::addConstantDensityFluid(
            builtin,request.composition,legacyFluidDiffusion);
    } else if (request.templateOrigin == PhysicsTemplateKind::EulerianEulerian) {
        Compose::addEulerianEulerianTemplate(builtin,result,request.phaseNames);
    } else if (request.pressureConstraint) {
        // Equation-source 请求：单流体压力约束方程族。它不是 solver family。
        Compose::addPressureConstraintFluid(builtin,*request.pressureConstraint);
    } else {
        // Legacy multiphase/level-set/thermal 模板：这些路径的守恒连续体来自
        // 同一个 authoritative preset。纯单流体必须显式选择 equation source。
        if (request.templateOrigin == PhysicsTemplateKind::SingleFluid) {
            throw std::runtime_error(
                "A single-fluid system must declare its equation source "
                "(single-fluid preset, declared equations, or a "
                "pressure-constraint request).");
        }
        Preset::installSingleFluid(
            builtin,SingleFluidPresetSpec{legacyFluidDiffusion});
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
    {
        SystemContribution contribution(result.rawSystem);
        Physics::SourceContribution::contribute(contribution,config.sources);
        models.applyContribution(std::move(contribution));
    }

    if (request.homogeneousThermodynamics) {
        models.recordContribution(
            "model.homogeneousMultiphase","homogeneous multiphase equations");
        addUnknown(models,"phaseMassAux","homogeneous phase mass auxiliary",1,
                   UnknownRole::Transported,
                   StorageBinding::SpecializedExecutor,{},0,"homogeneous");
        addEquation(models,{
            "E_PHASE_MASS","homogeneous phase-mass transport",
            "conservation",{"phaseMassAux"}},
            Equation::named("E_PHASE_MASS",
                Equation::ddt({"phaseMassAux"})
                    + Equation::div({"phaseMassFlux"})
                    == Equation::Symbol{"phaseMassSources"}));
    } else if (request.legacyMixture) {
        models.recordContribution(
            "model.legacyMixture","legacy mixture equations");
        addUnknown(models,"alphaAux","legacy transported volume fraction",1,
                   UnknownRole::Transported,
                   StorageBinding::SpecializedExecutor,{},0,"legacyMultiphase");
        addEquation(models,{
            "E_LEGACY_ALPHA","legacy volume-fraction transport",
            "conservation",{"alphaAux"}},
            Equation::named("E_LEGACY_ALPHA",
                Equation::ddt({"alphaAux"})
                    + Equation::div({"alphaFlux"})
                    == Equation::Symbol{"zero"}));
    }
    if (request.levelSet) {
        SystemContribution contribution(result.rawSystem);
        Physics::InterfaceModels::LevelSetContribution::contribute(
            contribution,*request.levelSet);
        models.applyContribution(std::move(contribution));
    }
    if (request.turbulence) {
        SystemContribution contribution(result.rawSystem);
        Turbulence::contribute(contribution,*request.turbulence);
        models.applyContribution(std::move(contribution));
    }
    if (request.immersed) {
        if (request.templateOrigin == PhysicsTemplateKind::EulerianEulerian) {
            if (request.immersed->enforcement
                == FDM::IBMEnforcement::GhostCell) {
                throw std::runtime_error(
                    "Eulerian multiphase + Ghost IBM is unsupported: "
                    "phase-wise ghost-state boundary closure is unavailable.");
            }
            throw std::runtime_error(
                "Eulerian multiphase IBM constraint exists, but required "
                "phase-wise IBM fluid-port assembly is unavailable.");
        }
        SystemContribution contribution(result.rawSystem);
        IBM::SystemContribution::contribute(
            contribution,*request.immersed);
        models.applyContribution(std::move(contribution));
    }
    // §15 precedence：用户修改必须排在 builtin defaults 与 model
    // contributions 之后、formulation/transformation 之前。当前没有 typed
    // lowering 的动作会显式失败，绝不静默覆盖。
    for (const SystemModification& modification : request.userModifications) {
        models.applyModification(modification);
    }

    TransformerRegistry transformers;
    transformers.registerTransformer(makePressureConstraintTransformer());
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
    if (request.immersed && request.immersed->monolithic
        && couplingPolicyContributed) {
        // DLM/KKT 是完全隐式约束耦合：它取代 segregated pressure policy，
        // 不是与之并存。preset 仍然必须被显式报告为未生效。
        result.executionPolicies.erase(
            std::remove_if(
                result.executionPolicies.begin(),result.executionPolicies.end(),
                [](const ExecutionPolicy& policy) {
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
    // STATE REALIZATION：只由 executable system 的 role/binding 导出。
    result.realization = compileStateRealization(result.executableSystem);
    // Time/constraint plan leaves are executable operations too. Their IDs are
    // declared before planning so provider resolution can validate every leaf.
    if (result.realization.conservativeTransportedMass) {
        const Provenance origin{OriginKind::Generated,"explicitTimeRecipe"};
        for (const char* id : {"flow.step.prepare", "flow.dt.compute",
                "flow.step.begin", "explicit.stage.execute",
                "flow.step.commit", "time.commit"}) {
            result.executableSystem.operations.push_back({
                id,id,OperationStage::Prepare,
                {OperationCapability::ConservativeExplicit},origin});
        }
    }
    for (const ExecutionPolicy& policy : result.executionPolicies) {
        const char* id = policy.kind == ExecutionPolicyKind::ConstraintProjection
            ? "ibm.constraint.project"
            : policy.kind == ExecutionPolicyKind::MonolithicKKT
                ? "ibm.kkt.solve" : nullptr;
        if (id) result.executableSystem.operations.push_back({
            id,id,OperationStage::Prepare,
            {OperationCapability::ImmersedConstraint},policy.origin});
    }
    // COUPLING PLAN CONTRIBUTION：片段只引用 formulation 声明的
    // executable operations；缺失的 stage 由 provider resolver 报告。
    std::vector<PlanFragment> planFragments;
    if (request.coupling && couplingPolicyContributed) {
        const auto sharedPolicy = std::find_if(
            result.executionPolicies.begin(),result.executionPolicies.end(),
            [](const ExecutionPolicy& policy) {
                return policy.id == kSharedPressureScheduleId;
            });
        PlanFragment fragment = sharedPolicy == result.executionPolicies.end()
            ? couplingPlanFragment(*request.coupling)
            : sharedPressurePlanFragment(*sharedPolicy,result.executableSystem);
        // monolithic KKT 已经显式取代该 preset（status=Inactive + reason）时
        // 保持那个报告，不把它降级成 provider 问题。
        const bool reportResolution =
            result.coupling.status == CouplingStatus::Active;
        std::vector<std::string> unresolved;
        if (!reportResolution) {
            // 保持既有状态报告（例如 superseded by monolithic KKT）。
        } else if (couplingPlanResolved(
                fragment,result.executableSystem,&unresolved)) {
            result.coupling.derivedOperations.clear();
            for (const ExecutableOperation& operation
                 : result.executableSystem.operations) {
                result.coupling.derivedOperations.push_back(
                    operation.operation);
            }
        } else {
            result.coupling.status = CouplingStatus::Unsupported;
            std::string list;
            for (const std::string& id : unresolved) {
                if (!list.empty()) list += ", ";
                list += id;
            }
            result.coupling.reason =
                "the registered coupling preset requires operations the "
                "formulation does not declare: [" + list + "]";
            result.coupling.derivedOperations = unresolved;
        }
        planFragments.push_back(std::move(fragment));
    }
    auto numericalRecipes = config.numerics.recipes;
    numericalRecipes.time = result.timeRecipe;
    const bool coreDiffusionConsumer = std::any_of(
        result.executableSystem.equations.begin(),
        result.executableSystem.equations.end(),
        [&](const EquationDescriptor& equation) {
            if (equation.id != "E_MASS" && equation.id != "E_MOMENTUM"
                && equation.id != "E_ENERGY") return false;
            const auto& definition = result.executableSystem
                .equationDefinitions.at(equation.id);
            const auto hasDiffusion = [](const Equation::Expression& expression) {
                return std::any_of(expression.terms.begin(),expression.terms.end(),
                    [](const Equation::Term& term) {
                        return term.kind == Equation::TermKind::Diffusion;
                    });
            };
            return hasDiffusion(definition.left)
                || hasDiffusion(definition.right);
        });
    if (coreDiffusionConsumer && !numericalRecipes.diffusion
        && !config.numerics.termRecipesDeclared) {
        numericalRecipes.diffusion =
            FDM::builtInDiffusionRecipe(config.numerics.viscous);
    }
    result.numericalSystem = NumericalCompiler::compile(
        result.executableSystem,numericalRecipes);
    // dt 限制策略属于 compiled HOW：runtime 读这份冻结值，而不是 raw config。
    result.numericalSystem.dt.cfl = config.numerics.cfl;
    result.numericalSystem.dt.maxDeltaT = config.numerics.maxDeltaT;
    result.numericalSystem.phaseTransport.convection =
        config.pressure.phaseTransport.convection;
    result.numericalSystem.phaseTransport.sourceCfl =
        config.pressure.phaseTransport.sourceCfl;
    for (const auto& provider : result.numericalSystem.providerRequirements) {
        requireProvider(result,provider,
                        "execute a compiled mathematical term binding");
    }
    result.solvePlan = SolvePlanner::compile(
        result.executableSystem,result.executionPolicies,result.timeRecipe,
        planFragments);
    result.runtime.capabilities = compileExecutionCapabilities(
        result.executableSystem,result.executionPolicies);
    result.runtime.operationBindings = resolveOperationBindings(
        result.executableSystem,result.realization,result.numericalSystem,
        result.solvePlan,result.executionPolicies);
    result.runtime.report = reportOperationBindings(
        result.solvePlan,result.runtime.operationBindings);

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
    const bool constraintDof = request.immersed
        && request.immersed->introducesMultiplier && request.parallel;
    result.runtime.requirements.push_back({
        "ConstraintGlobalDof",constraintDof,
        !constraintDof || request.capabilities.canonicalConstraintDof,
        "unique owner for Lambda/lambda and sparse J/S edges"});
    const bool distributedSolve = request.immersed
        && request.immersed->monolithic;
    result.runtime.requirements.push_back({
        "DistributedLinearSystem",distributedSolve,
        !distributedSolve || request.capabilities.distributedLinearSystem,
        "GlobalDofId is mapped to backend rows outside equation assembly"});
    const bool homogeneousUnsupported = request.homogeneousThermodynamics
        && (!result.realization.conservativeTransportedMass
            || (config.turbulence.enabled
                && config.turbulence.family != FDM::TurbulenceFamily::DNS));
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
    const bool unsupportedPhaseChange = request.phaseChange
        && (request.levelSet || request.legacyMixture);
    result.runtime.requirements.push_back({
        "PhaseChangeExecution",request.phaseChange,
        !unsupportedPhaseChange,
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
    return std::any_of(system.unknowns.begin(),system.unknowns.end(),
        [&](const UnknownDescriptor& value) { return value.id == id; });
}

bool hasUnknown(const ExecutableEquationSystem& system, std::string_view id) {
    return std::any_of(system.unknowns.begin(),system.unknowns.end(),
        [&](const UnknownDescriptor& value) { return value.id == id; });
}

bool hasEquation(const ResolvedSimulationSystem& system, std::string_view id) {
    return hasEquation(system.executableSystem,id);
}

bool hasEquation(const RawEquationSystem& system, std::string_view id) {
    return std::any_of(system.equations.begin(),system.equations.end(),
        [&](const EquationDescriptor& value) { return value.id == id; });
}

bool hasEquation(const ExecutableEquationSystem& system, std::string_view id) {
    return std::any_of(system.equations.begin(),system.equations.end(),
        [&](const EquationDescriptor& value) { return value.id == id; });
}

const Equation::Definition& equationDefinition(
        const ResolvedSimulationSystem& system, std::string_view id) {
    return system.executableSystem.equationDefinitions.at(std::string(id));
}

bool hasEquationPrefix(
        const ExecutableEquationSystem& system, std::string_view prefix) {
    return std::any_of(
        system.equations.begin(),
        system.equations.end(),
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
