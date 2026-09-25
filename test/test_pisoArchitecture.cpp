#include "solver/system/SF_systemBuilder.h"
#include "solver/system/SF_systemValidator.h"
#include "solver/system/SF_solvePlan.h"
#include "solver/run/SF_planExecutor.h"
#include "core/interfaces/SF_immersedSystem.h"
#include "SF_pressureCoupling.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

using namespace SF;

[[noreturn]] void fail(const std::string& message) {
    std::cerr << "PISO architecture test failed: " << message << '\n';
    std::exit(1);
}

void require(bool condition, const std::string& message) {
    if (!condition) fail(message);
}

FDM::SolverConfig pisoConfig() {
    FDM::SolverConfig config;
    config.numerics.timeRecipe =
        FDM::builtInTimeRecipe(FDM::TimeRecipeId::ForwardEuler);
    config.pressure.coupling.preset = FDM::PressureCouplingPreset::PISO;
    config.pressure.coupling.pressureCorrectors = 1;
    config.pressure.coupling.nonOrthogonalCorrectors = 0;
    return config;
}

const System::CompiledEquation& compiled(
        const System::ResolvedSimulationSystem& system,
        const std::string& id) {
    const auto found = std::find_if(
        system.executableSystem.compiledEquations.begin(),
        system.executableSystem.compiledEquations.end(),
        [&](const System::CompiledEquation& item) {
            return item.equationId == id;
        });
    if (found == system.executableSystem.compiledEquations.end()) {
        fail("missing compiled equation " + id);
    }
    return *found;
}

void collectOperations(
        const System::SolvePlanNode& node,
        std::vector<System::OpId>& result) {
    if (!node.operation.empty()) {
        result.push_back(node.operation);
    }
    for (const auto& child : node.children) collectOperations(child,result);
}

void collectStageLoops(
        const System::SolvePlanNode& node,
        std::vector<const System::SolvePlanNode*>& result) {
    if (node.kind == System::PlanNodeKind::StageLoop) result.push_back(&node);
    for (const auto& child : node.children) collectStageLoops(child,result);
}

const System::SolvePlanNode* findNode(
        const System::SolvePlanNode& node,
        const std::string& id) {
    if (node.id == id) return &node;
    for (const auto& child : node.children) {
        if (const auto* found = findNode(child,id)) return found;
    }
    return nullptr;
}

void requireOperationLeaves(const System::SolvePlanNode& node) {
    if (node.children.empty()) {
        require(!node.operation.empty(),
                "compiled Plan contains a leaf without an OpId: "+node.id);
    }
    for (const auto& child : node.children) requireOperationLeaves(child);
}

std::string rawSignature(const System::RawEquationSystem& raw) {
    std::ostringstream out;
    for (const auto& unknown : raw.unknowns) {
        out << "U:" << unknown.id << ':' << unknown.components << ':'
            << static_cast<int>(unknown.role) << ':' << unknown.storageKey << '\n';
    }
    for (const auto& equation : raw.equations) {
        out << "E:" << equation.id << ':' << static_cast<int>(equation.category);
        for (const auto& unknown : equation.solvedUnknowns) out << ':' << unknown;
        out << '\n';
        const auto& definition = raw.equationDefinitions.at(equation.id);
        for (const auto* expression : {&definition.left,&definition.right}) {
            for (const auto& term : expression->terms) {
                out << "T:" << static_cast<int>(term.kind) << ':'
                    << term.primary.name << ':' << term.secondary.name << '\n';
            }
            out << "|\n";
        }
    }
    for (const auto& constraint : raw.constraints) {
        out << "C:" << constraint.id << ':' << constraint.equation << ':'
            << constraint.multiplierUnknown << '\n';
    }
    for (const auto& closure : raw.closures) out << "K:" << closure << '\n';
    return out.str();
}

} // namespace

int main() {
    auto providers = builtinCompositionProviders();
    providers.registerEquationOfState("dummyEos");
    providers.registerAlgorithm("dummyAlgorithm");
    require(providers.hasEquationOfState("dummyEos")
                && providers.hasAlgorithm("dummyAlgorithm"),
            "open composition registry requires a central enum edit");

    // coupling preset 由 case 输入注册；它必须与同一份 workflow 同步（生产
    // 路径中二者都来自 config.solver.pressure）。
    const auto requestFor = [](const FDM::SolverConfig& config) {
        System::BuildRequest item;
        item.templateOrigin = System::PhysicsTemplateKind::SingleFluid;
        item.pressureConstraint = System::PressureConstraintSpec{false};
        item.coupling = System::couplingRequestFrom(
            config.pressure.coupling,true);
        return item;
    };
    System::BuildRequest request = requestFor(pisoConfig());
    const auto system = System::build(pisoConfig(),request);
    System::validate(system);

    require(!System::hasUnknown(system.rawSystem,"pPrime"),
            "pPrime leaked into RawEquationSystem");
    require(!System::hasEquation(system.rawSystem,"E_PRESSURE"),
            "pressure-correction equation leaked into RawEquationSystem");
    require(System::hasEquation(system.rawSystem,"E_MOMENTUM"),
            "raw momentum equation is absent");
    require(System::hasConstraint(system.rawSystem,"C_INCOMPRESSIBILITY"),
            "raw incompressibility constraint is absent");

    require(System::hasUnknown(system,"pPrime"),
            "transformer did not generate pPrime workspace");
    require(System::hasEquation(system,"E_MOMENTUM_PREDICTOR"),
            "transformer did not generate momentum predictor");
    require(System::hasEquation(system,"E_PRESSURE"),
            "transformer did not generate pressure correction");

    const auto& predictor = compiled(system,"E_MOMENTUM_PREDICTOR");
    require(predictor.operatorBinding
                == "momentum.predictor",
            "momentum predictor operator binding is wrong");
    require(predictor.assemblesRhs && !predictor.assemblesMatrix,
            "momentum predictor matrix/RHS contract is wrong");
    const auto& pressure = compiled(system,"E_PRESSURE");
    require(pressure.operatorBinding
                == "pressure.correction",
            "pressure operator binding is wrong");
    require(pressure.assemblesMatrix && pressure.assemblesRhs,
            "pressure matrix/RHS contract is incomplete");
    require(std::any_of(
                pressure.resources.begin(),pressure.resources.end(),
                [](const System::CompiledResourceBinding& resource) {
                    return resource.storage == "pressureMatrix"
                        && resource.access == System::ResourceAccessMode::Write;
                }),
            "pressure matrix storage binding is absent");

    require(system.runtime.report.status == System::RuntimeStatus::Runnable,
            "minimal PISO was not compiled as a runnable plan");
    require(system.solvePlan.root.kind == System::PlanNodeKind::Sequence,
            "PISO root is not a sequence");
    const auto* outer = findNode(
        system.solvePlan.root,"PISO.outerCorrectors");
    const auto* loop = findNode(
        system.solvePlan.root,"PISO.pressureCorrectors");
    const auto* nonOrthogonal = findNode(
        system.solvePlan.root,"PISO.nonOrthogonalCorrectors");
    require(outer && loop && nonOrthogonal,
            "pressure-corrector loop is absent");
    require(outer->repetitions == 1 && loop->repetitions == 1
                && nonOrthogonal->repetitions == 1,
            "pressure-corrector loop count is not compiled from policy");
    requireOperationLeaves(system.solvePlan.root);

    std::vector<System::OpId> operations;
    collectOperations(system.solvePlan.root,operations);
    const std::vector<System::OpId> expected = {
        "pressure.prepare", "momentum.assemble", "momentum.solve",
        "pressure.boundary.prepare", "pressure.assemble", "pressure.solve",
        "velocity.correct", "pressure.update.prepare", "flux.correct",
        "pressure.correction.commit", "pressure.step.commit"};
    require(operations == expected,
            "structured PISO operation order differs from compiled plan contract");
    require(system.runtime.report.requiredOperations == expected,
            "runtime operation requirements are not recursively Plan-derived");
    require(System::requiresProvider(system,"flow.conservative")
                && !System::requiresProvider(system,"flow.eulerian-pressure"),
            "single-fluid PISO selected a solver-family execution route");
    require(System::requiresProvider(system,"linear.hypre")
                && System::requiresRuntimeService(system,"mpi.initialized"),
            "PISO provider requirements omit its HYPRE/MPI runtime dependency");

    System::CompiledSolvePlan genericPlan;
    genericPlan.root = {System::PlanNodeKind::Sequence,"root","root",{}, {},
                        {},1,{
        {System::PlanNodeKind::Update,"a","a",{}, {},"op.a",1,{}},
        {System::PlanNodeKind::Loop,"loop","loop",{}, {},{},3,{
            {System::PlanNodeKind::Update,"b","b",{}, {},"op.b",1,{}}}}}};
    std::vector<std::string> trace;
    Run::OpRegistry genericOps;
    genericOps.bind("op.a",[&] { trace.push_back("a"); });
    genericOps.bind("op.b",[&] { trace.push_back("b"); });
    Run::PlanExecutor::execute(genericPlan,genericOps);
    require(trace == std::vector<std::string>{"a","b","b","b"},
            "open plan executor did not preserve sequence/loop order");

    System::CompiledSolvePlan stagePlan;
    stagePlan.root = {System::PlanNodeKind::StageLoop,"stages","stages",
                      {}, {},{},4,{
        {System::PlanNodeKind::Update,"stage","stage",{}, {},
         "test.stage",1,{}}}};
    std::vector<std::pair<int,int>> stages;
    Run::OpRegistry stageOps;
    stageOps.bind("test.stage",[&](const Run::ExecutionContext& context) {
        stages.push_back({context.stageIndex,context.stageCount});
    });
    Run::PlanExecutor::execute(stagePlan,stageOps);
    require(stages == std::vector<std::pair<int,int>>{
                {0,4},{1,4},{2,4},{3,4}},
            "StageLoop did not expose the complete zero-based stage context");
    require(System::SolvePlanner::requiredOperations(genericPlan)
                == std::vector<System::OpId>{"op.a","op.b"},
            "recursive operation discovery lost a nested Plan leaf");
    System::CompiledSolvePlan missingPlan;
    missingPlan.root = {System::PlanNodeKind::Sequence,"missing","missing",
                        {}, {},{},1,{
        {System::PlanNodeKind::Update,"prepare","prepare",{}, {},
         "flow.step.prepare",1,{}},
        {System::PlanNodeKind::Update,"ibm","ibm",{}, {},
         "ibm.constraint.project",1,{}}}};
    bool missingIdVisible = false;
    try {
        Run::OpRegistry incomplete;
        incomplete.bind("flow.step.prepare",[] {});
        Run::PlanExecutor::validateBindings(missingPlan,incomplete);
    } catch (const std::runtime_error& error) {
        missingIdVisible = std::string(error.what()).find(
            "ibm.constraint.project") != std::string::npos;
    }
    require(missingIdVisible,
            "missing operation provider did not fail before execution");
    Run::OpRegistry complete;
    complete.bind("flow.step.prepare",[] {});
    complete.bind("ibm.constraint.project",[] {});
    Run::PlanExecutor::validateBindings(missingPlan,complete);
    Run::PlanExecutor::execute(missingPlan,complete);

    std::string referenceRawSystem;
    for (const auto& item : std::vector<std::pair<FDM::TimeRecipeId,int>>{
            {FDM::TimeRecipeId::ForwardEuler,1},
            {FDM::TimeRecipeId::SSPRK3,3},
            {FDM::TimeRecipeId::ClassicalRK4,4}}) {
        FDM::SolverConfig explicitConfig;
        explicitConfig.numerics.timeRecipe = FDM::builtInTimeRecipe(item.first);
        System::BuildRequest explicitRequest;
        explicitRequest.templateOrigin =
            System::PhysicsTemplateKind::SingleFluid;
        explicitRequest.singleFluidPreset =
            System::SingleFluidPresetSpec{false};
        const auto explicitSystem =
            System::build(explicitConfig,explicitRequest);
        const std::string signature = rawSignature(explicitSystem.rawSystem);
        if (referenceRawSystem.empty()) referenceRawSystem = signature;
        require(signature == referenceRawSystem,
                "RawEquationSystem changed with the selected TimeRecipe");
        require(explicitSystem.timeRecipe.id() == item.first
                    && explicitSystem.timeRecipe.stageCount() == item.second,
                "resolved system lost the immutable TimeRecipe contract");
        std::vector<const System::SolvePlanNode*> stageLoops;
        collectStageLoops(explicitSystem.solvePlan.root,stageLoops);
        require(stageLoops.size() == 1,
                "explicit policies did not aggregate into one StageLoop");
        require(stageLoops.front()->repetitions == item.second,
                "compiled explicit StageLoop has the wrong repetition count");
        require(System::SolvePlanner::requiredOperations(explicitSystem.solvePlan)
                    == std::vector<System::OpId>{
                        "flow.step.prepare","flow.dt.compute","flow.step.begin",
                        "explicit.stage.execute","flow.step.commit","time.commit"},
                "explicit plan does not expose the truthful fused operation order");
    }

    auto twoCorrectorConfig = pisoConfig();
    twoCorrectorConfig.pressure.coupling.pressureCorrectors = 2;
    const auto twoCorrector =
        System::build(twoCorrectorConfig,requestFor(twoCorrectorConfig));
    require(twoCorrector.runtime.report.status == System::RuntimeStatus::Runnable,
            "PISO corrector count was not routed by the capability signature");
    require(twoCorrector.runtime.capabilities.pressureCorrectors == 2,
            "compiled capability signature lost the PISO corrector count");
    require(twoCorrector.runtime.capabilities.outerCorrectors == 1
                && twoCorrector.runtime.capabilities.nonOrthogonalCorrectors == 0,
            "pressure schedule count meanings are inconsistent");

    auto renamedPolicies = system.executionPolicies;
    for (auto& policy : renamedPolicies) {
        if (policy.id == "S_PRESSURE") policy.strategyName = "display-only";
    }
    // Planner 只做 merge/order：coupling preset 的 plan fragment 提供顺序与
    // 重复次数，stage -> OpId 由 executable operation authority 解析。
    const auto renamedPlan = System::SolvePlanner::compile(
        system.executableSystem,renamedPolicies,system.timeRecipe,
        {System::couplingPlanFragment(*request.coupling)});
    require(System::SolvePlanner::requiredOperations(renamedPlan)
                == expected,
            "strategyName still controls pressure Plan lowering");

    for (const auto algorithm : {
            FDM::PressureCouplingPreset::SIMPLE,
            FDM::PressureCouplingPreset::PIMPLE}) {
        auto fixedPointConfig = pisoConfig();
        fixedPointConfig.pressure.coupling.preset = algorithm;
        fixedPointConfig.pressure.coupling.outerCorrectors = 2;
        fixedPointConfig.pressure.coupling.pressureCorrectors = 3;
        fixedPointConfig.pressure.coupling.nonOrthogonalCorrectors = 1;
        const auto fixedPoint =
            System::build(fixedPointConfig,requestFor(fixedPointConfig));
        System::validate(fixedPoint);
        const auto policy = std::find_if(
            fixedPoint.executionPolicies.begin(),
            fixedPoint.executionPolicies.end(),
            [](const System::ExecutionPolicy& item) {
                return item.id == "S_PRESSURE";
            });
        require(policy != fixedPoint.executionPolicies.end()
                    && policy->kind
                        == System::ExecutionPolicyKind::PressureVelocityFixedPoint
                    && policy->strategyKind
                        == FDM::SolveStrategyKind::PressureVelocityCoupling,
                "SIMPLE/PIMPLE did not retain typed fixed-point semantics");
        require(policy->repeatCount == 2
                    && policy->nestedRepeatCount == 3
                    && policy->innerRepeatCount == 2,
                "single-fluid pressure policy lost outer/pressure/non-orthogonal counts");
        require(fixedPoint.runtime.report.status
                    == System::RuntimeStatus::Unsupported
                    && !fixedPoint.runtime.report.missingOperations.empty(),
                "fixed-point pressure schedule was reported Runnable");
        require(fixedPoint.runtime.report.reason.find("full dt")
                    != std::string::npos,
                "fixed-point provider failure does not explain the predictor guard");
        const auto fixedOps = System::SolvePlanner::requiredOperations(
            fixedPoint.solvePlan);
        require(std::find(fixedOps.begin(),fixedOps.end(),"momentum.solve")
                    == fixedOps.end()
                    && std::find(fixedOps.begin(),fixedOps.end(),
                                 "pressure.schedule.predictor.solve")
                        != fixedOps.end(),
                "fixed-point schedule reuses the physical-time momentum solve");
        const auto* fixedOuter = findNode(
            fixedPoint.solvePlan.root,"PressureSchedule.outerCorrectors");
        const auto* fixedPressure = findNode(
            fixedPoint.solvePlan.root,"PressureSchedule.pressureCorrectors");
        const auto* fixedNonOrthogonal = findNode(
            fixedPoint.solvePlan.root,
            "PressureSchedule.nonOrthogonalCorrectors");
        require(fixedOuter && fixedPressure && fixedNonOrthogonal
                    && fixedOuter->repetitions == 2
                    && fixedPressure->repetitions == 3
                    && fixedNonOrthogonal->repetitions == 2,
                "fixed-point Plan does not preserve all three schedule levels");
        requireOperationLeaves(fixedPoint.solvePlan.root);
    }

    auto eulerianConfig = pisoConfig();
    eulerianConfig.pressure.coupling.preset = FDM::PressureCouplingPreset::PIMPLE;
    eulerianConfig.pressure.coupling.outerCorrectors = 2;
    eulerianConfig.pressure.coupling.pressureCorrectors = 3;
    eulerianConfig.pressure.coupling.nonOrthogonalCorrectors = 1;
    System::BuildRequest eulerianRequest;
    eulerianRequest.templateOrigin =
        System::PhysicsTemplateKind::EulerianEulerian;
    eulerianRequest.phaseNames = {"water","air"};
    // Eulerian 的共享压力 schedule 也由注册的 coupling preset 贡献。
    eulerianRequest.coupling = System::couplingRequestFrom(
        eulerianConfig.pressure.coupling,true);
    const auto eulerian = System::build(eulerianConfig,eulerianRequest);
    require(System::requiresProvider(eulerian,"flow.eulerian-pressure")
                && !System::requiresProvider(eulerian,"flow.conservative"),
            "Eulerian operation group did not select its phase-state provider");
    const auto eulerianOps =
        System::SolvePlanner::requiredOperations(eulerian.solvePlan);
    const std::vector<System::OpId> frozenEulerianOps{
        "ee.dt.compute", "ee.step.begin", "ee.interphase.compute",
        "ee.sources.assemble", "ee.sources.validate",
        "ee.momentum.diagonal", "ee.momentum.flux",
        "ee.faceFlux.canonical", "ee.continuity.assemble",
        "ee.boundary.prepare", "ee.momentum.solve",
        "ee.interphase.correct", "ee.boundary.afterMomentum",
        "ee.diagonal.sync", "ee.momentum.flux.after",
        "ee.faceFlux.canonical.after", "ee.pressure.solve",
        "ee.pressure.publish", "ee.pressure.sync", "ee.phase.correct",
        "ee.faceFlux.correct", "ee.boundary.afterPressure",
        "ee.faceFlux.canonical.pressure", "ee.energy.solve",
        "ee.boundary.final", "ee.outer.validate", "ee.step.commit",
        "ee.time.commit"};
    require(eulerianOps == frozenEulerianOps,
            "generic Eulerian fragment lowering changed operation order");
    const auto* eulerianOuter = findNode(eulerian.solvePlan.root,"EE.outer");
    const auto* eulerianPressure = findNode(eulerian.solvePlan.root,"EE.pressure");
    const auto* eulerianNonOrthogonal = findNode(
        eulerian.solvePlan.root,"EE.nonOrthogonal");
    require(eulerian.solvePlan.root.id == "EE.step"
                && eulerianOuter && eulerianOuter->repetitions == 2
                && eulerianPressure && eulerianPressure->repetitions == 3
                && eulerianNonOrthogonal
                && eulerianNonOrthogonal->repetitions == 2,
            "generic Eulerian fragment lowering changed nested loop topology");
    require(std::find(eulerianOps.begin(),eulerianOps.end(),
                      "ee.turbulence.solve") == eulerianOps.end(),
            "Eulerian Plan scheduled turbulence without turbulence equations");
    bool unconsumedVisible = false;
    try {
        auto policies = eulerian.executionPolicies;
        System::ExecutionPolicy orphan;
        orphan.id = "S_ORPHAN";
        orphan.name = "unconsumed test contribution";
        orphan.kind = System::ExecutionPolicyKind::BoundaryClosure;
        orphan.strategyName = "test";
        orphan.strategyKind = FDM::SolveStrategyKind::BoundaryClosure;
        policies.push_back(orphan);
        (void)System::SolvePlanner::compile(
            eulerian.executableSystem,policies,eulerian.timeRecipe);
    } catch (const std::runtime_error& error) {
        unconsumedVisible = std::string(error.what()).find("S_ORPHAN")
            != std::string::npos;
    }
    require(unconsumedVisible,
            "Eulerian planning silently ignored an active policy");
    FDM::ImmersedAlgorithmDescriptor eulerianImmersed;
    eulerianImmersed.id = "testEulerianIBM";
    eulerianImmersed.enforcement = FDM::IBMEnforcement::FractionalDLM;
    auto unsupportedEulerianRequest = eulerianRequest;
    unsupportedEulerianRequest.immersed = &eulerianImmersed;
    bool eulerianIbmRejected = false;
    try {
        (void)System::build(eulerianConfig,unsupportedEulerianRequest);
    } catch (const std::runtime_error& error) {
        eulerianIbmRejected = std::string(error.what()).find(
            "phase-wise IBM fluid-port assembly is unavailable")
            != std::string::npos;
    }
    require(eulerianIbmRejected,
            "Eulerian variational IBM was silently omitted from execution");

    auto threePhaseRequest = eulerianRequest;
    threePhaseRequest.phaseNames = {"water","air","oil"};
    const auto threePhase = System::build(eulerianConfig,threePhaseRequest);
    require(System::requiresProvider(threePhase,"flow.eulerian-pressure")
                && System::requiresProvider(threePhase,"linear.hypre"),
            "N-phase state changed execution-provider group selection");

    auto invalidProviderCoverage = system;
    invalidProviderCoverage.runtime.providerRequirements.clear();
    bool missingProviderVisible = false;
    try {
        System::validate(invalidProviderCoverage);
    } catch (const std::runtime_error& error) {
        missingProviderVisible = std::string(error.what()).find(
            "Resolved operation provider is absent from runtime requirements")
            != std::string::npos;
    }
    require(missingProviderVisible,
            "missing flow provider requirement did not fail during validation");

    System::BuildRequest constantDensityRequest;
    constantDensityRequest.templateOrigin = System::PhysicsTemplateKind::SingleFluid;
    constantDensityRequest.composition.declared = true;
    constantDensityRequest.composition.equations = {"Momentum","Continuity"};
    constantDensityRequest.composition.thermoDynamics.equationOfState =
        "rhoConst";
    constantDensityRequest.composition.thermoDynamics.transport =
        "const";
    constantDensityRequest.composition.thermoDynamics.constantDensity = 1.0;
    constantDensityRequest.composition.algorithm = "PISO";
    constantDensityRequest.composition.pressureCorrectors = 2;
    constantDensityRequest.coupling = System::couplingRequestFrom(
        pisoConfig().pressure.coupling,true);
    const auto constantDensity = System::build(pisoConfig(),constantDensityRequest);
    System::validate(constantDensity);
    require(System::hasEquation(constantDensity.rawSystem,"E_CONTINUITY"),
            "rhoConst composition did not keep Continuity as a raw equation");
    require(!System::hasUnknown(constantDensity.rawSystem,"rhoE"),
            "rhoConst composition fabricated a total-energy unknown");
    require(constantDensity.classification.densityBehavior == "constant"
                && constantDensity.classification.thermodynamicCompressibility
                       == "zero",
            "rhoConst closure was not resolved as constant density");
    require(constantDensity.runtime.report.status == System::RuntimeStatus::Unsupported,
            "rhoConst composition silently selected a legacy pressure runtime");
    const auto bindingFor = [](const System::ResolvedSimulationSystem& resolved,
                               const System::OpId& id)
            -> const System::ResolvedOperationBinding* {
        const auto& bindings = resolved.runtime.operationBindings;
        const auto found = std::find_if(bindings.begin(),bindings.end(),
            [&](const System::ResolvedOperationBinding& binding) {
                return binding.operation == id;
            });
        return found == bindings.end() ? nullptr : &*found;
    };
    const auto* conservativePressure = bindingFor(system,"pressure.solve");
    const auto* constantPressure = bindingFor(constantDensity,"pressure.solve");
    require(conservativePressure && constantPressure
                && conservativePressure->status == System::BindingStatus::Resolved
                && conservativePressure->provider == "flow.conservative"
                && constantPressure->status == System::BindingStatus::Unsupported
                && constantPressure->provider.empty()
                && constantPressure->reason.find("constant-density pressure-multiplier")
                    != std::string::npos,
            "pressure.solve did not resolve by capability and state realization");
    std::vector<System::OpId> unresolved;
    for (const auto& binding : constantDensity.runtime.operationBindings) {
        if (binding.status == System::BindingStatus::Unsupported) {
            unresolved.push_back(binding.operation);
        }
    }
    require(unresolved == constantDensity.runtime.report.missingOperations,
            "RuntimeReport missing operations differ from unresolved bindings");

    System::BuildRequest sodRequest;
    sodRequest.templateOrigin = System::PhysicsTemplateKind::SingleFluid;
    sodRequest.singleFluidPreset = System::SingleFluidPresetSpec{false};
    sodRequest.composition.declared = true;
    sodRequest.composition.equations = {"Continuity","Momentum","Energy"};
    sodRequest.composition.thermoDynamics.equationOfState =
        "perfectGas";
    sodRequest.composition.thermoDynamics.thermo = "hConst";
    sodRequest.composition.algorithm = "Explicit";
    const auto sod = System::build(pisoConfig(),sodRequest);
    System::validate(sod);
    require(sod.classification.densityBehavior == "variable"
                && sod.classification.thermodynamicCompressibility
                       == "available",
            "perfectGas composition did not resolve variable-density closure");
    require(System::hasUnknown(sod.rawSystem,"rhoE")
                && !System::hasConstraint(sod.rawSystem,"C_INCOMPRESSIBILITY"),
            "Sod composition was polluted by a pressure constraint");
    require(System::requiresProvider(sod,"flow.conservative")
                && System::requiresProvider(
                    sod,"thermodynamics.single-fluid"),
            "explicit single-fluid system lacks derived numerical providers");

    auto levelSetRequest = sodRequest;
    levelSetRequest.levelSet =
        Physics::InterfaceModels::LevelSetContribution::Spec{};
    const auto levelSet = System::build(pisoConfig(),levelSetRequest);
    require(System::requiresProvider(levelSet,"equation.level-set"),
            "level-set equation did not derive its state/equation provider");

    std::cout << "PISO transformation/compiled-binding/plan contract passed\n";
    return 0;
}
