#include "solver/system/SF_systemBuilder.h"
#include "solver/system/SF_systemValidator.h"
#include "solver/system/SF_solvePlan.h"
#include "solver/system/SF_termProviderCatalog.h"
#include "solver/run/SF_planExecutor.h"
#include "core/interfaces/SF_immersedSystem.h"
#include "models/physics/interfaceModel/levelSet/SF_levelSetSystemContribution.h"
#include "models/physics/SF_sourceContribution.h"
#include "SF_pressureCoupling.h"
#include "solver/system/SF_methodObjects.h"

#include <algorithm>
#include <cmath>
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

struct TestBodyForceModel {
    double strength;

    System::SystemContribution contribution() const {
        System::SystemContribution value;
        value.recordContribution("model.testBodyForce","test body force");
        value.extendMathematics({"momentum"},System::FormulaExpr::op("source",
            {System::FormulaExpr::symbol("customForce")},"customForce"));
        return value;
    }

    System::SourceTermProviderDescriptor provider() const {
        return {"source.customForce.primitive",
            [](const System::TermMatchContext& context) {
                return context.expression.kind==System::FormulaExpr::Kind::Operator
                    && context.expression.name=="source"
                    && context.expression.arguments.size()==1
                    && context.expression.arguments.front().name=="customForce"
                    && context.output=="U";
            },{},[value=strength] {
                return System::PrimitiveMomentumSource{
                    [value](const Field&,int,int,int,const Vector3&) {
                        return Vector3(value,0.0,0.0);
                    }};
            },"test.bodyForce"};
    }
};

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
    for (const auto& unknown : raw.state.symbols()) {
        out << "U:" << unknown.id << ':' << unknown.components << ':'
            << static_cast<int>(unknown.role) << ':' << unknown.storageKey << '\n';
    }
    for (const auto& equation : raw.legacyEquations) {
        out << "E:" << equation.id << ':' << static_cast<int>(equation.category);
        for (const auto& unknown : equation.solvedUnknowns) out << ':' << unknown;
        out << '\n';
        const auto& definition = raw.legacyDefinitions.at(equation.id);
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
    item.composition.stateDeclared = true;
    item.composition.solutionVariables = {"rho","rhoU","rhoE"};
        item.templateOrigin = System::PhysicsTemplateKind::SingleFluid;
        item.pressureConstraint = System::PressureConstraintSpec{false};
        item.coupling = System::couplingRequestFrom(
            config.pressure.coupling,true);
        return item;
    };
    System::BuildRequest request = requestFor(pisoConfig());
    const auto system = System::build(pisoConfig(),request);
    System::validate(system);
    auto invalidPiso=*request.coupling;
    invalidPiso.outerCorrectors=2;
    const auto invalidPisoMatch=System::matchPressureCoupling(
        invalidPiso,system.rawSystem);
    require(invalidPisoMatch.status==System::CouplingStatus::Invalid
                && invalidPisoMatch.reason.find("use PIMPLE")!=std::string::npos,
            "PISO silently accepted an outer fixed-point iteration");
    auto invalidSimple=*request.coupling;
    invalidSimple.presetKind=FDM::PressureCouplingPreset::SIMPLE;
    invalidSimple.preset="SIMPLE";
    invalidSimple.pressureCorrectors=2;
    const auto invalidSimpleMatch=System::matchPressureCoupling(
        invalidSimple,system.rawSystem);
    require(invalidSimpleMatch.status==System::CouplingStatus::Invalid
                && invalidSimpleMatch.reason.find("use PIMPLE")!=std::string::npos,
            "SIMPLE silently accepted multiple pressure corrections");

    require(!System::hasUnknown(system.rawSystem,"pPrime"),
            "pPrime leaked into RawEquationSystem");
    require(!System::hasEquation(system.rawSystem,"pSimple"),
            "pressure-correction equation leaked into RawEquationSystem");
    require(System::hasEquation(system.rawSystem,"momentum"),
            "raw momentum equation is absent");
    require(System::hasConstraint(system.rawSystem,"C_INCOMPRESSIBILITY"),
            "raw incompressibility constraint is absent");

    require(!System::hasUnknown(system,"pPrime") && !System::hasUnknown(system,"phi"),
            "conservative pressure invented a base correction/face-flux authority");
    require(system.executableSystem.state.at("rhoU").storageKey=="conservative"
                && system.executableSystem.state.at("U").derivation==System::StateDerivation::Velocity,
            "conservative momentum was mislabeled as primitive velocity");
    require(system.executionPolicies.empty() && system.executableSystem.compiledEquations.empty(),
            "single-fluid pressure retained legacy schedule/equation authority");
    const auto& calls=system.solvePlan.compiledProgram.steps;
    require(calls.size()==6 && calls.front().equationMethod=="ConservativePressureMomentum"
                && calls.front().sourceMathInputs==std::vector<std::string>{"continuity","energy"}
                && calls.front().writes==std::vector<std::string>{"rho","rhoU","rhoE"}
                && calls.front().calls.size()==3,
            "native conservative predictor lacks an explicit mass/momentum/energy fusion contract");
    require(calls[1].target.symbol=="p" && calls[1].target.kind==System::TargetKind::Correction
                && calls[3].target.symbol=="p" && calls[3].target.kind==System::TargetKind::Working
                && calls[4].target.kind==System::TargetKind::Workspace
                && calls[5].target.symbol=="rhoE",
            "native conservative pressure targets do not match correction/EOS publication mathematics");
    for (const auto& call:calls) require(call.backendProvider=="flow.conservative",
            "native conservative method did not freeze its backend owner");
    const auto owned=[&](const auto& self,const System::SolvePlanNode& node)->void {
        if (!node.operation.empty()) require(!node.legacyAdapter && node.provider=="flow.conservative",
            "single-fluid pressure leaf still depends on a legacy provider selector");
        for (const auto& child:node.children) self(self,child);
    };
    owned(owned,system.solvePlan.root);
    for (const auto& view:system.solvePlan.compiledProgram.stateViews) {
        if (view.kind==System::StateViewKind::Correction || view.kind==System::StateViewKind::Working)
            require(view.symbol=="p" && view.owner==System::StateViewOwner::NumericalProvider,
                "pressure views acquired a second compiler workspace authority");
    }

    require(system.runtime.report.status == System::RuntimeStatus::Runnable,
            "minimal PISO was not compiled as a runnable plan");
    require(system.solvePlan.root.kind == System::PlanNodeKind::Sequence,
            "PISO root is not a sequence");
    const auto* outer = findNode(
        system.solvePlan.root,"outer");
    const auto* loop = findNode(
        system.solvePlan.root,"pressure");
    const auto* nonOrthogonal = findNode(
        system.solvePlan.root,"nonOrthogonal");
    require(!outer && loop && nonOrthogonal,
            "pressure-corrector loop is absent");
    require(loop->repetitions == 1 && nonOrthogonal->repetitions == 1,
            "pressure-corrector loop count is not compiled from policy");
    requireOperationLeaves(system.solvePlan.root);

    std::vector<System::OpId> operations;
    collectOperations(system.solvePlan.root,operations);
    const std::vector<System::OpId> expected = {
        "pressure.prepare", "momentum.assemble", "momentum.solve",
        "pressure.boundary.prepare", "pressure.assemble", "pressure.solve",
        "velocity.correct", "pressure.update.prepare", "flux.correct",
        "pressure.correction.commit", "pressure.step.commit", "time.commit"};
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

    // The executor knows only a signal ID. A callback decides convergence;
    // the complete iteration body runs before the generic Loop may exit.
    System::CompiledSolvePlan terminatingPlan;
    auto& root=terminatingPlan.root;
    root.kind=System::PlanNodeKind::Sequence;
    root.id="mock.step";
    System::SolvePlanNode begin;
    begin.kind=System::PlanNodeKind::Update;
    begin.id="mock.begin";
    begin.operation="mock.begin";
    root.children.push_back(begin);
    System::SolvePlanNode iterations;
    iterations.kind=System::PlanNodeKind::Loop;
    iterations.id="mock.outer";
    iterations.repetitions=5;
    iterations.terminationSignal="mock.converged";
    for (const auto& id:{"mock.predict","mock.evaluate","mock.end"}) {
        System::SolvePlanNode leaf;
        leaf.kind=System::PlanNodeKind::Update;
        leaf.id=id;
        leaf.operation=id;
        iterations.children.push_back(std::move(leaf));
    }
    root.children.push_back(std::move(iterations));
    System::SolvePlanNode commit;
    commit.kind=System::PlanNodeKind::Commit;
    commit.id="mock.commit";
    commit.operation="mock.commit";
    root.children.push_back(commit);
    int iteration=0,target=2,commits=0;
    double physicalTime=0.0,base=2.0,working=base,observedBase=base;
    Run::OpRegistry terminatingOps;
    terminatingOps.bind("mock.begin",[&](const Run::ExecutionContext& context) {
        iteration=0;
        observedBase=base;
        working=base;
        context.signals->reset("mock.converged");
    });
    terminatingOps.bind("mock.predict",[&] {
        ++iteration;
        working=observedBase+0.5*(4.0-working);
    });
    terminatingOps.bind("mock.evaluate",[&](const Run::ExecutionContext& context) {
        context.signals->publish("mock.converged",iteration>=target);
    });
    int ended=0;
    terminatingOps.bind("mock.end",[&] { ++ended; });
    terminatingOps.bind("mock.commit",[&] {
        ++commits;
        base=working;
        physicalTime+=0.5;
    });
    Run::PlanExecutor::execute(terminatingPlan,terminatingOps);
    require(iteration==2 && ended==2 && commits==1
                && std::abs(working-2.5)<1e-14
                && std::abs(physicalTime-0.5)<1e-14,
            "generic early exit changed the fixed physical-time base or commit");
    target=3;
    Run::PlanExecutor::execute(terminatingPlan,terminatingOps);
    require(iteration==3 && ended==5 && commits==2
                && std::abs(physicalTime-1.0)<1e-14,
            "loop termination signal leaked into the next physical step");

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
    explicitRequest.composition.stateDeclared = true;
    explicitRequest.composition.solutionVariables = {"rho","rhoU","rhoE"};
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

    // A different relation cannot be accepted by a fixed existing kernel.
    auto changedMath=system.executableSystem;
    auto changed=changedMath.registry.at("correctP");
    changed.rhs=System::FormulaExpr::constantValue(1.0);
    changedMath.registry.replace(changed);
    System::ExecutionProgram authoredHow;
    authoredHow.root=system.solvePlan.compiledProgram.root;
    bool changedMathRejected=false;
    try {
        (void)System::compileExecutionProgram(changedMath,authoredHow,
            system.numericalSelection.bindings,System::builtinProviders(),
            &system.numericalSystem.time.recipe);
    } catch (const std::runtime_error& error) {
        changedMathRejected=std::string(error.what()).find("changed mathematics")!=std::string::npos;
    }
    require(changedMathRejected,"conservative pressure provider accepted unsupported authored mathematics");

    for (const auto algorithm : {
            FDM::PressureCouplingPreset::SIMPLE,
            FDM::PressureCouplingPreset::PIMPLE}) {
        auto fixedPointConfig = pisoConfig();
        fixedPointConfig.pressure.coupling.preset = algorithm;
        fixedPointConfig.pressure.coupling.outerCorrectors = 2;
        fixedPointConfig.pressure.coupling.pressureCorrectors =
            algorithm==FDM::PressureCouplingPreset::SIMPLE ? 1 : 3;
        fixedPointConfig.pressure.coupling.nonOrthogonalCorrectors = 1;
        const auto fixedPoint =
            System::build(fixedPointConfig,requestFor(fixedPointConfig));
        System::validate(fixedPoint);
        require(fixedPoint.executionPolicies.empty(),
                "SIMPLE/PIMPLE conservative HOW retained a duplicate legacy schedule");
        require(fixedPoint.runtime.report.status
                    == System::RuntimeStatus::Unsupported
                    && !fixedPoint.runtime.report.missingOperations.empty(),
                "fixed-point pressure schedule was reported Runnable");
        require(fixedPoint.runtime.report.reason.find(
                    "pressure.step.begin") != std::string::npos,
                "conservative fixed-point path did not report its missing "
                "fixed-time provider");
        const auto fixedOps = System::SolvePlanner::requiredOperations(
            fixedPoint.solvePlan);
        std::vector<System::OpId> fixedOrder;
        collectOperations(fixedPoint.solvePlan.root,fixedOrder);
        const std::vector<System::OpId> expectedFixedOrder = {
            "pressure.prepare","pressure.step.begin","pressure.iteration.begin",
            "momentum.assemble","momentum.solve",
            "pressure.boundary.prepare","pressure.assemble","pressure.solve",
            "velocity.correct","pressure.update.prepare","flux.correct",
            "pressure.correction.commit","pressure.relaxation.apply",
            "pressure.flux.consistency.restore","pressure.convergence.evaluate",
            "pressure.iteration.end","pressure.step.commit","time.commit"};
        require(fixedOrder==expectedFixedOrder,
            "SIMPLE/PIMPLE flattened pressure operation order changed");
        require(std::find(fixedOps.begin(),fixedOps.end(),"momentum.solve")
                    != fixedOps.end()
                    && std::find(fixedOps.begin(),fixedOps.end(),
                                 "pressure.step.begin")!=fixedOps.end()
                    && std::find(fixedOps.begin(),fixedOps.end(),
                                 "pressure.relaxation.apply")!=fixedOps.end()
                    && std::find(fixedOps.begin(),fixedOps.end(),
                                 "pressure.flux.consistency.restore")
                        !=fixedOps.end(),
                "fixed-point schedule lost its frozen-time operation contract");
        const auto* fixedOuter = findNode(
            fixedPoint.solvePlan.root,"outer");
        const auto* fixedPressure = findNode(
            fixedPoint.solvePlan.root,"pressure");
        const auto* singleCorrection = findNode(
            fixedPoint.solvePlan.root,"pressureCorrection");
        const auto* fixedNonOrthogonal = findNode(
            fixedPoint.solvePlan.root,
            "nonOrthogonal");
        require(fixedOuter && fixedNonOrthogonal
                    && fixedOuter->repetitions == 2
                    && fixedOuter->terminationSignal
                        == System::kPressureOuterConvergedSignal
                    && (algorithm==FDM::PressureCouplingPreset::SIMPLE
                        ? singleCorrection && !fixedPressure
                        : fixedPressure && fixedPressure->repetitions==3)
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
    eulerianRequest.composition.stateDeclared = true;
    eulerianRequest.composition.solutionVariables = {"rho","rhoU","rhoE"};
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
        System::LegacyExecutionPolicy orphan;
        orphan.id = "S_ORPHAN";
        orphan.name = "unconsumed test contribution";
        orphan.kind = System::LegacyExecutionPolicyKind::BoundaryClosure;
        orphan.strategyName = "test";
        orphan.strategyKind = FDM::SolveStrategyKind::BoundaryClosure;
        policies.push_back(orphan);
        (void)System::SolvePlanner::compile(
            eulerian.executableSystem,policies,
            eulerian.numericalSystem.time.recipe);
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
    unsupportedEulerianRequest.unsupportedCompositionReason =
        "Eulerian multiphase IBM constraint exists, but required "
        "phase-wise IBM fluid-port assembly is unavailable.";
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
    constantDensityRequest.composition.stateDeclared = true;
    constantDensityRequest.composition.solutionVariables = {"rho","rhoU","rhoE"};
    constantDensityRequest.templateOrigin = System::PhysicsTemplateKind::SingleFluid;
    constantDensityRequest.composition.declared = true;
    constantDensityRequest.composition.solutionVariables = {"U","p"};
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
    require(!System::hasEquation(constantDensity.rawSystem,"E_CONTINUITY")
                && System::hasConstraint(constantDensity.rawSystem,
                                         "C_INCOMPRESSIBILITY"),
            "rhoConst composition must lower Continuity to div(U)=0, not transport rho");
    require(!constantDensity.rawSystem.state.contains("rhoE")
                && !System::hasEquation(constantDensity.rawSystem,"energy"),
            "rhoConst composition activated a catalog total-energy state/equation");
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

    // A complete primitive numerical recipe resolves the same mathematical
    // PISO plan. WENO remains unsupported for U/p; it must not silently select
    // primitive upwind or the legacy conservative corrector.
    auto primitiveConfig=pisoConfig();
    primitiveConfig.numerics.formulation=
        FDM::EquationFormulation::PrimitiveDifferential;
    primitiveConfig.numerics.convection=FDM::ConvectionScheme::Upwind1;
    primitiveConfig.numerics.reconstruction=FDM::ReconstructionVariable::Primitive;
    primitiveConfig.numerics.flux=FDM::FluxSplitter::UpwindAdvection;
    primitiveConfig.numerics.recipes.convection=
        FDM::TermRecipe::primitiveUpwindRecipe();
    const auto primitive=System::build(primitiveConfig,constantDensityRequest);
    System::validate(primitive);
    require(primitive.runtime.report.status==System::RuntimeStatus::Runnable,
            "primitiveUpwind1 constant-density PISO did not resolve all operations");
    require(bindingFor(primitive,"pressure.solve")->provider
                == "flow.pressure-operators"
                && bindingFor(primitive,"flux.correct")->provider
                    == "flow.rhie-chow",
            "constant-density PISO did not bind neutral pressure/face-flux operations");
    auto heatOnlyRequest=constantDensityRequest;
    FDM::SourceConfig heatConfig;
    heatConfig.enabled={FDM::SourceKind::WallHeat};
    System::SystemContribution heatContribution;
    Physics::SourceContribution::contribute(heatContribution,heatConfig);
    heatOnlyRequest.modelContributions.push_back(std::move(heatContribution));
    bool heatRejected=false;
    try {
        (void)System::build(primitiveConfig,heatOnlyRequest);
    } catch (const std::runtime_error& error) {
        heatRejected=std::string(error.what()).find(
            "Source contribution has no registered target equation")!=std::string::npos;
    }
    require(heatRejected,
            "wallHeat on pressure-only equations was silently ignored");
    auto customRequest=constantDensityRequest;
    {
        const TestBodyForceModel model{0.01};
        customRequest.termProviders.push_back(model.provider());
        customRequest.modelContributions.push_back(model.contribution());
    }
    const auto customSystem=System::build(primitiveConfig,customRequest);
    System::validate(customSystem);
    require(customSystem.runtime.report.status==System::RuntimeStatus::Runnable,
            "typed custom Momentum term with registered provider is not runnable");
    require(std::any_of(customSystem.numericalSystem.operators.begin(),
                        customSystem.numericalSystem.operators.end(),
        [](const System::CompiledSpatialBinding& term) {
            return term.primary=="customForce"
                && term.executionEquationId=="momentum"
                && term.provider=="source.customForce.primitive"
                && term.providerOwner=="test.bodyForce"
                && term.compiledDataAvailable
                && term.primitiveSource;
        }),"custom term did not bind to the transformed Momentum predictor");
    const auto compiledForce=std::find_if(
        customSystem.numericalSystem.operators.begin(),
        customSystem.numericalSystem.operators.end(),
        [](const auto& term) { return term.primary=="customForce"; });
    Field testGeometry;
    const auto executed=compiledForce->primitiveSource(
        testGeometry,0,0,0,Vector3());
    require(executed.x==0.01 && executed.y==0.0 && executed.z==0.0,
            "test model did not execute after its configuration lifetime");
    std::vector<System::OpId> baseOperations,customOperations;
    collectOperations(primitive.solvePlan.root,baseOperations);
    collectOperations(customSystem.solvePlan.root,customOperations);
    require(baseOperations==customOperations,
            "custom Momentum source changed pressure solve-plan operations");
    for (const auto recipe:{FDM::TimeRecipeId::SSPRK3,
                           FDM::TimeRecipeId::ClassicalRK4}) {
        for (const auto preset:{FDM::PressureCouplingPreset::PISO,
                               FDM::PressureCouplingPreset::SIMPLE,
                               FDM::PressureCouplingPreset::PIMPLE}) {
            auto unsupportedConfig=primitiveConfig;
            unsupportedConfig.numerics.timeRecipe=FDM::builtInTimeRecipe(recipe);
            auto unsupportedRequest=constantDensityRequest;
            unsupportedRequest.composition.algorithm=FDM::toString(preset);
            unsupportedRequest.coupling->presetKind=preset;
            unsupportedRequest.coupling->preset=FDM::toString(preset);
            bool rejected=false;
            try {
                const auto unsupported=System::build(unsupportedConfig,unsupportedRequest);
                rejected=unsupported.runtime.report.status!=System::RuntimeStatus::Runnable;
            } catch (const std::runtime_error& error) {
                rejected=std::string(error.what()).find("Unsupported: provider")!=std::string::npos;
            }
            require(rejected,"multi-stage pressure capability was accepted or silently replaced");
        }
    }
    for (const auto preset:{FDM::PressureCouplingPreset::SIMPLE,
                            FDM::PressureCouplingPreset::PIMPLE}) {
        auto fixed=constantDensityRequest;
        fixed.composition.algorithm=FDM::toString(preset);
        fixed.coupling->presetKind=preset;
        fixed.coupling->preset=FDM::toString(preset);
        fixed.coupling->outerCorrectors=2;
        const auto selected=System::build(primitiveConfig,fixed);
        require(selected.runtime.report.status==System::RuntimeStatus::Runnable,
                "constant-density fixed-time pressure schedule is not runnable");
        for (const System::OpId& operation:{
                "pressure.assemble", "pressure.solve",
                "velocity.correct", "flux.correct"}) {
            const auto* pisoBinding=bindingFor(primitive,operation);
            const auto* fixedBinding=bindingFor(selected,operation);
            require(pisoBinding && fixedBinding
                        && pisoBinding->status==System::BindingStatus::Resolved
                        && fixedBinding->status==System::BindingStatus::Resolved
                        && fixedBinding->provider==pisoBinding->provider,
                    "fixed-point plan selected a different provider for "
                        +operation);
        }
        const auto operations=System::SolvePlanner::requiredOperations(
            selected.solvePlan);
        require(std::find(operations.begin(),operations.end(),
                          "pressure.step.begin")!=operations.end()
                    && std::find(operations.begin(),operations.end(),
                                 "pressure.relaxation.apply")!=operations.end(),
                "fixed-point plan lacks base capture or relaxation");
    }
    auto distributedPressure = constantDensityRequest;
    distributedPressure.parallel = true;
    const auto distributed = System::build(primitiveConfig,distributedPressure);
    require(distributed.runtime.report.status == System::RuntimeStatus::Runnable
                && distributed.runtime.report.missingOperations.empty(),
            "distributed constant-density pressure lost its shared operations");
    auto nonOrthogonalPressure = constantDensityRequest;
    nonOrthogonalPressure.coupling->nonOrthogonalCorrectors = 1;
    const auto unsupportedNonOrthogonal = System::build(
        primitiveConfig,nonOrthogonalPressure);
    require(unsupportedNonOrthogonal.runtime.report.status
                == System::RuntimeStatus::Unsupported
                && !unsupportedNonOrthogonal.runtime.report.missingOperations.empty(),
            "non-orthogonal constant-density correction was reported Runnable");

    System::BuildRequest sodRequest;
    sodRequest.composition.stateDeclared = true;
    sodRequest.composition.solutionVariables = {"rho","rhoU","rhoE"};
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

    require(!sod.executableSystem.state.contains("k") && !sod.executableSystem.state.contains("omega")
            && !sod.executableSystem.state.contains("alpha") && !sod.executableSystem.state.contains("h")
            && !primitive.executableSystem.state.contains("rhoE"),
            "Unused catalog symbols leaked into active case STATE.");
    auto requestedState=sodRequest;
    System::SystemContribution closureRequest;
    closureRequest.requireState("h");
    requestedState.modelContributions.push_back(closureRequest);
    const auto withEnthalpy=System::build(pisoConfig(),requestedState);
    require(withEnthalpy.executableSystem.state.contains("h")
            && withEnthalpy.executableSystem.registry.entries().size()==sod.executableSystem.registry.entries().size()
            && withEnthalpy.solvePlan.compiledProgram.steps.size()==sod.solvePlan.compiledProgram.steps.size()
            && !withEnthalpy.executableSystem.state.contains("k"),
            "Module STATE request required repeated metadata or activated equations/models.");
    bool unknownRequestRejected=false;
    closureRequest.requiredStates={"customWithoutMetadata"};
    requestedState.modelContributions={closureRequest};
    try { (void)System::build(pisoConfig(),requestedState); }
    catch (const std::runtime_error& error) {
        unknownRequestRejected=std::string(error.what()).find("addState")!=std::string::npos;
    }
    require(unknownRequestRejected,"Unknown custom STATE request inferred metadata.");

    auto levelSetRequest = sodRequest;
    levelSetRequest.modelRequiresDiffusion=true;
    levelSetRequest.additionalExecutionContributions=true;
    System::SystemContribution levelSetContribution;
    Physics::InterfaceModels::LevelSetContribution::contribute(
        levelSetContribution,{});
    levelSetRequest.modelContributions.push_back(
        std::move(levelSetContribution));
    const auto levelSet = System::build(pisoConfig(),levelSetRequest);
    require(System::requiresProvider(levelSet,"equation.level-set"),
            "level-set equation did not derive its state/equation provider");

    std::cout << "PISO transformation/compiled-binding/plan contract passed\n";
    return 0;
}
