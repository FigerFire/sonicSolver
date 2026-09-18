#include "solver/system/SF_systemBuilder.h"
#include "solver/system/SF_systemValidator.h"
#include "solver/system/SF_solvePlan.h"
#include "solver/run/SF_planExecutor.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>
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
    config.numerics.solver = FDM::SolverAlgorithm::PressureBased;
    config.numerics.time = FDM::TimeScheme::Euler;
    config.pressure.workflow.type = FDM::SolverAlgorithm::PressureBased;
    config.pressure.workflow.algorithm = FDM::PressureAlgorithm::PISO;
    config.pressure.workflow.pressureCorrectors = 1;
    config.pressure.workflow.nonOrthogonalCorrectors = 0;
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

} // namespace

int main() {
    auto providers = builtinCompositionProviders();
    providers.registerEquationOfState("dummyEos");
    providers.registerAlgorithm("dummyAlgorithm");
    require(providers.hasEquationOfState("dummyEos")
                && providers.hasAlgorithm("dummyAlgorithm"),
            "open composition registry requires a central enum edit");

    System::BuildRequest request;
    request.templateOrigin = System::PhysicsTemplateKind::SingleFluid;
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

    require(system.runtime.status == System::RuntimeStatus::Runnable,
            "minimal PISO was not compiled as a runnable plan");
    require(system.runtime.requiredAdapters.empty(),
            "minimal PISO still reports a legacy execution adapter");
    require(system.solvePlan.root.kind == System::PlanNodeKind::Sequence,
            "PISO root is not a sequence");
    const auto loop = std::find_if(
        system.solvePlan.root.children.begin(),system.solvePlan.root.children.end(),
        [](const System::SolvePlanNode& node) {
            return node.kind == System::PlanNodeKind::Loop
                && node.id == "PISO.pressureCorrectors";
        });
    require(loop != system.solvePlan.root.children.end(),
            "pressure-corrector loop is absent");
    require(loop->repetitions == 1,
            "pressure-corrector loop count is not compiled from policy");

    std::vector<System::OpId> operations;
    collectOperations(system.solvePlan.root,operations);
    const std::vector<System::OpId> expected = {
        "pressure.prepare", "momentum.assemble", "momentum.solve",
        "pressure.boundary.prepare", "pressure.assemble", "pressure.solve",
        "velocity.correct", "pressure.update.prepare", "flux.correct",
        "pressure.correction.commit", "pressure.step.commit"};
    require(operations == expected,
            "structured PISO operation order differs from compiled plan contract");
    require(system.runtime.requiredOperations == expected,
            "runtime operation requirements are not recursively Plan-derived");

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
    bool missingIdVisible = false;
    try {
        System::CompiledSolvePlan missingPlan;
        missingPlan.root = {System::PlanNodeKind::Update,"missing","missing",
                            {}, {},"missing.operation",1,{}};
        Run::PlanExecutor::execute(missingPlan,genericOps);
    } catch (const std::runtime_error& error) {
        missingIdVisible = std::string(error.what()).find("missing.operation")
            != std::string::npos;
    }
    require(missingIdVisible,"missing operation ID did not fail visibly");

    for (const auto& item : std::vector<std::pair<FDM::TimeScheme,int>>{
            {FDM::TimeScheme::Euler,1},
            {FDM::TimeScheme::SSPRK3,3},
            {FDM::TimeScheme::RK4,4}}) {
        FDM::SolverConfig explicitConfig;
        explicitConfig.numerics.solver = FDM::SolverAlgorithm::DensityBased;
        explicitConfig.numerics.time = item.first;
        System::BuildRequest explicitRequest;
        explicitRequest.templateOrigin =
            System::PhysicsTemplateKind::SingleFluid;
        const auto explicitSystem =
            System::build(explicitConfig,explicitRequest);
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
    twoCorrectorConfig.pressure.workflow.pressureCorrectors = 2;
    const auto twoCorrector = System::build(twoCorrectorConfig,request);
    require(twoCorrector.runtime.status == System::RuntimeStatus::Runnable,
            "PISO corrector count was not routed by the capability signature");
    require(twoCorrector.executionCapabilities.pressureCorrectors == 2,
            "compiled capability signature lost the PISO corrector count");

    auto eulerianConfig = pisoConfig();
    eulerianConfig.pressure.workflow.algorithm = FDM::PressureAlgorithm::PIMPLE;
    eulerianConfig.pressure.workflow.outerCorrectors = 2;
    eulerianConfig.pressure.workflow.pressureCorrectors = 3;
    eulerianConfig.pressure.workflow.nonOrthogonalCorrectors = 1;
    System::BuildRequest eulerianRequest;
    eulerianRequest.templateOrigin =
        System::PhysicsTemplateKind::EulerianEulerian;
    eulerianRequest.phaseNames = {"water","air"};
    const auto eulerian = System::build(eulerianConfig,eulerianRequest);
    const auto eulerianOps =
        System::SolvePlanner::requiredOperations(eulerian.solvePlan);
    require(!eulerianOps.empty()
                && eulerianOps.front() == "ee.dt.compute"
                && eulerianOps.back() == "ee.time.commit",
            "Eulerian Plan does not own dt calculation and clock commit");
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
            eulerian.executableSystem,policies,eulerian.timeIntegrator);
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
    const auto constantDensity = System::build(pisoConfig(),constantDensityRequest);
    System::validate(constantDensity);
    require(System::hasEquation(constantDensity.rawSystem,"E_CONTINUITY"),
            "rhoConst composition did not keep Continuity as a raw equation");
    require(!System::hasUnknown(constantDensity.rawSystem,"rhoE"),
            "rhoConst composition fabricated a total-energy unknown");
    require(constantDensity.densityBehavior == "constant"
                && constantDensity.thermodynamicCompressibility == "zero",
            "rhoConst closure was not resolved as constant density");
    require(constantDensity.runtime.status == System::RuntimeStatus::Unsupported,
            "rhoConst composition silently selected a legacy pressure runtime");

    System::BuildRequest sodRequest;
    sodRequest.templateOrigin = System::PhysicsTemplateKind::SingleFluid;
    sodRequest.composition.declared = true;
    sodRequest.composition.equations = {"Continuity","Momentum","Energy"};
    sodRequest.composition.thermoDynamics.equationOfState =
        "perfectGas";
    sodRequest.composition.thermoDynamics.thermo = "hConst";
    sodRequest.composition.algorithm = "Explicit";
    const auto sod = System::build(pisoConfig(),sodRequest);
    System::validate(sod);
    require(sod.densityBehavior == "variable"
                && sod.thermodynamicCompressibility == "available",
            "perfectGas composition did not resolve variable-density closure");
    require(System::hasUnknown(sod.rawSystem,"rhoE")
                && !System::hasConstraint(sod.rawSystem,"C_INCOMPRESSIBILITY"),
            "Sod composition was polluted by a pressure constraint");

    std::cout << "PISO transformation/compiled-binding/plan contract passed\n";
    return 0;
}
