#include "solver/system/SF_methodObjects.h"
#include "solver/run/SF_planExecutor.h"
#include "solver/system/SF_pressureCoupling.h"
#include "solver/system/SF_solvePlan.h"
#include <algorithm>
#include "core/system/SF_operationIds.h"

#include <stdexcept>
#include <string>

namespace {
void require(bool value,const char* message) {
    if (!value) throw std::runtime_error(message);
}

SF::System::StateSymbol stored(std::string id) {
    SF::System::StateSymbol result;
    result.id=std::move(id);
    result.storageKey="testStorage";
    return result;
}

class SyntheticGroup final:public SF::System::IProvider {
public:
    std::string name="SyntheticGroup",owner="test.group";
    std::vector<SF::System::CompiledMathRef> members;
    std::string_view id() const override { return name; }
    std::string_view runtimeProvider() const override { return owner; }
    SF::System::CompiledEquationCall compile(const SF::System::ExecutableEquationSystem&,
            const SF::System::EquationCall& call,const SF::System::NumericalBinding&) const override {
        using namespace SF::System;
        CompiledEquationCall result;
        result.source=call;result.target.symbol=call.target.symbol;result.target.kind=call.target.kind;
        result.equationMethod=name;result.calls={{call.equation,call.target.symbol}};
        result.fusionKey="explicitGroup";result.fusionMembers=members;
        for (const auto* operation:{"test.prepare","test.evaluate"}) {
            SolvePlanNode leaf;leaf.kind=PlanNodeKind::Update;leaf.id=operation;
            leaf.operation=operation;leaf.equationCalls=result.calls;
            result.fragment.children.push_back(std::move(leaf));
        }
        return result;
    }
};

void checkGenericFusion() {
    using namespace SF::System;
    for (const int size:{2,3}) {
        ExecutableEquationSystem system;
        SyntheticGroup method;
        ExecutionProgram program;
        std::vector<NumericalBinding> bindings;
        // Deliberately reverse registry insertion: the explicit contract owns order.
        for (int i=size-1;i>=0;--i) {
            const auto state="value"+std::to_string(i),equation="relation"+std::to_string(i);
            system.state.add(stored(state));
            system.registry.add({equation,FormulaExpr::symbol(state),FormulaExpr::constantValue(0),{}});
        }
        system.state.add(stored("otherValue"));
        for (int i=0;i<size;++i) {
            const auto state="value"+std::to_string(i),equation="relation"+std::to_string(i);
            method.members.push_back({equation,state});
            ExecutionScope call;call.kind=ExecutionKind::EquationCall;call.step={equation,{state}};
            program.root.children.push_back(call);
            bindings.push_back({equation,std::string(method.id())});
        }
        SyntheticGroup other=method;other.name="OtherGroup";other.owner="test.other";
        ProviderRegistry providers;providers.add(method);providers.add(other);
        const auto compiled=compileExecutionProgram(system,program,bindings,providers);
        require(compiled.steps.size()==static_cast<std::size_t>(size)
            && compiled.loweredRoot.children.size()==1,
            "Generic explicit group lost source calls or duplicated a numerical fragment.");
        for (const auto& leaf:compiled.loweredRoot.children.front().children) {
            require(leaf.equationCalls.size()==static_cast<std::size_t>(size)
                && leaf.provider=="test.group","Generic fusion lost member provenance or owner.");
            for (int i=0;i<size;++i)
                require(leaf.equationCalls[i].equation==method.members[i].equation
                    && leaf.equationCalls[i].target==method.members[i].target,
                    "Registry insertion reordered an explicit numerical group.");
        }
        const auto reject=[&](ExecutionProgram candidate,std::vector<NumericalBinding> selected) {
            bool rejected=false;
            try { (void)compileExecutionProgram(system,candidate,selected,providers); }
            catch (const std::runtime_error&) { rejected=true; }
            require(rejected,"Generic compiler accepted an invalid explicit fusion group.");
        };
        auto invalid=program;invalid.root.children.pop_back();
        auto shortened=bindings;shortened.pop_back();reject(invalid,shortened);
        invalid=program;invalid.root.children.back().step.equation=method.members.front().equation;
        shortened=bindings;shortened.pop_back();reject(invalid,shortened);
        invalid=program;invalid.root.children.back().step.target.symbol="otherValue";reject(invalid,bindings);
        invalid=program;std::swap(invalid.root.children.front(),invalid.root.children.back());reject(invalid,bindings);
        invalid=program;ExecutionScope nested;nested.children={invalid.root.children.back()};
        invalid.root.children.back()=nested;reject(invalid,bindings);
        auto mixed=bindings;mixed.back().method=other.name;reject(program,mixed);
        const auto last=method.members.back();method.members.back()=method.members.front();
        invalid=program;invalid.root.children.back().step=invalid.root.children.front().step;
        shortened=bindings;shortened.pop_back();reject(invalid,shortened);
        method.members.back()=last;
    }
}
}

int main() {
    checkGenericFusion();
    using namespace SF::System;
    using Expr=FormulaExpr;
    ExecutableEquationSystem system;
    system.state.add(stored("C"));
    system.registry.add({"MyScalar",
        Expr::add(Expr::op("ddt",{Expr::symbol("C")}),
                  Expr::op("div",{Expr::symbol("scalarFlux")})),
        Expr::constantValue(0.0),{}});
    system.registry.add({"UpdateC",Expr::symbol("C"),
                         Expr::add(Expr::symbol("C"),Expr::constantValue(1.0)),{}});
    system.registry.add({"DiffuseC",
        Expr::negate(Expr::op("div",{
            Expr::multiply(Expr::symbol("nu"),
                Expr::op("grad",{Expr::symbol("C")}))},"diffusion")),
        Expr::symbol("S"),{}});
    ExecutionProgram program{EquationCall{"MyScalar",{"C"}}};
    const std::vector<NumericalBinding> binding{
        {"MyScalar","ConservativeResidual"}};
    const auto methods=builtinProviders();
    const auto compiled=compileExecutionProgram(system,program,binding,methods);
    require(compiled.steps.size()==1
            && compiled.steps.front().calls.size()==1
            && compiled.steps.front().calls.front().equation=="MyScalar"
            && compiled.steps.front().calls.front().target=="C"
            && compiled.steps.front().equationMethod=="ConservativeResidual",
            "Arbitrary FormulaGroup did not lower from explicit HOW output.");
    const auto temporalMethods=builtinTemporalMethods();
    const auto euler=temporalMethods.at(SF::FDM::TimeRecipeId::ForwardEuler).compile(
        SF::FDM::builtInTimeRecipe(SF::FDM::TimeRecipeId::ForwardEuler));
    const auto rk4=temporalMethods.at(SF::FDM::TimeRecipeId::ClassicalRK4).compile(
        SF::FDM::builtInTimeRecipe(SF::FDM::TimeRecipeId::ClassicalRK4));
    require(euler.stageCount()==1 && rk4.stageCount()==4
            && program.root.children.front().step.equation=="MyScalar"
            && compiled.steps.front().calls.front().target=="C",
            "TemporalMethod selection changed WHAT or HOW.");
    const ExecutionProgram assignment{EquationCall{"UpdateC",{"C"}}};
    const auto direct=compileExecutionProgram(system,assignment,
        {{"UpdateC","DirectEvaluation"}},methods);
    require(direct.steps.front().equationMethod=="DirectEvaluation",
            "DirectEvaluation did not lower its own numerical mode.");
    const auto linear=compileExecutionProgram(system,
        ExecutionProgram{EquationCall{"DiffuseC",{"C"}}},
        {{"DiffuseC","LinearEquation"}},methods);
    require(linear.steps.front().equationMethod=="LinearEquation"
            && linear.steps.front().requirements.size()==3,
            "LinearEquation did not supply its assembly/backend contract.");
    StateRegistry invalidState;
    bool qualifiedRegistrationRejected=false;
    try { invalidState.add(stored("U*")); }
    catch (const std::runtime_error&) { qualifiedRegistrationRejected=true; }
    require(qualifiedRegistrationRejected,"STATE accepted a HOW-qualified base symbol.");
    const auto& rkMethod=temporalMethods.at(SF::FDM::TimeRecipeId::ClassicalRK4);
    const auto rk=rkMethod.compile(SF::FDM::builtInTimeRecipe(SF::FDM::TimeRecipeId::ClassicalRK4));
    const auto staged=compileExecutionProgram(system,program,binding,methods,&rk,&rkMethod);
    require(std::count_if(staged.stateViews.begin(),staged.stateViews.end(),[](const auto& view) {
                return view.kind==StateViewKind::Stage;
            })==4 && staged.stateViews.size()==6 && system.state.size()==1,
            "WHICH RK4 did not lazily request stage views from the same base STATE.");

    // Pressure methods lower five explicit HOW outputs to validated numerical
    // fragments. Only the existing callbacks execute gauge/boundary/face math.
    ExecutableEquationSystem pressureSystem;
    for (const auto* id:{"U","p","pPrime","phi"})
        pressureSystem.state.add(stored(id));
    pressureSystem.legacyEquations.push_back(
        {"momentum","momentum","physical",{"U"}});
    pressureSystem.registry.add({"momentum",Expr::symbol("U"),
                                 Expr::symbol("U"),{}});
    pressureSystem.registry.add({"pSimple",
        Expr::op("grad",{Expr::symbol("pPrime")}),Expr::constantValue(0),{}});
    pressureSystem.registry.add({"correctU",Expr::symbol("U"),
                                 Expr::symbol("pPrime"),{}});
    pressureSystem.registry.add({"correctP",Expr::symbol("p"),
                                 Expr::symbol("pPrime"),{}});
    pressureSystem.registry.add({"correctFluxp",Expr::symbol("phi"),
                                 Expr::symbol("pPrime"),{}});
    ExecutionProgram pressureProgram{{
        {"momentum",{"U"}},
        {"pSimple",{"p",TargetKind::Correction}},
        {"correctU",{"U"}},
        {"correctP",{"p"}},
        {"correctFluxp",{"phi",TargetKind::Physical}}}};
    const std::vector<NumericalBinding> pressureBindings{
        {"momentum","PressureMomentum"},
        {"pSimple","PressureCorrection"},
        {"correctU","VelocityCorrection"},
        {"correctP","PressureUpdate"},
        {"correctFluxp","FluxCorrection"}};
    const auto pressureCompiled=compileExecutionProgram(
        pressureSystem,pressureProgram,pressureBindings,methods);
    const auto fragmentOperations=[](const CompiledEquationCall& step) {
        std::vector<OpId> result;
        for (const auto& child:step.fragment.children)
            result.push_back(child.operation);
        return result;
    };
    require(pressureCompiled.steps.size()==5
        && fragmentOperations(pressureCompiled.steps.at(0))
            ==std::vector<OpId>{OpIds::MomentumAssemble,OpIds::MomentumSolve}
        && fragmentOperations(pressureCompiled.steps.at(1))
            ==std::vector<OpId>{OpIds::PressureBoundaryPrepare,
                OpIds::PressureAssemble,OpIds::PressureSolve}
        && fragmentOperations(pressureCompiled.steps.at(2))
            ==std::vector<OpId>{OpIds::VelocityCorrect}
        && fragmentOperations(pressureCompiled.steps.at(3))
            ==std::vector<OpId>{OpIds::PressureUpdatePrepare}
        && fragmentOperations(pressureCompiled.steps.at(4))
            ==std::vector<OpId>{OpIds::FluxCorrect,OpIds::PressureCorrectionCommit},
        "Pressure EquationMethods lost their validated numerical fragments.");
    require(pressureCompiled.steps.at(1).target.workspace
                =="pressureCorrection"
            && pressureCompiled.steps.at(1).equationMethod
                =="PressureCorrection"
            && pressureCompiled.steps.at(4).requirements.size()==3,
        "Pressure method output, implicit solve, or canonical-face contract changed.");
    CompiledSolvePlan methodPlan;
    methodPlan.root=compileMethodProgram(pressureCompiled);
    std::vector<std::string> pressureOrder;
    SF::Run::OpRegistry pressureOperations;
    for (const auto* id:{OpIds::MomentumAssemble,OpIds::MomentumSolve,
            OpIds::PressureBoundaryPrepare,OpIds::PressureAssemble,
            OpIds::PressureSolve,OpIds::VelocityCorrect,
            OpIds::PressureUpdatePrepare,OpIds::FluxCorrect,OpIds::PressureCorrectionCommit,OpIds::PressurePrepare})
        pressureOperations.bind(id,id==OpIds::FluxCorrect ? "flow.rhie-chow" : "flow.pressure-operators",[&pressureOrder,id] { pressureOrder.push_back(id); });
    SF::Run::PlanExecutor::execute(methodPlan,pressureOperations);
    require(pressureOrder==std::vector<std::string>{
        OpIds::PressurePrepare,OpIds::MomentumAssemble,OpIds::MomentumSolve,
        OpIds::PressureBoundaryPrepare,OpIds::PressureAssemble,
        OpIds::PressureSolve,OpIds::VelocityCorrect,
        OpIds::PressureUpdatePrepare,OpIds::FluxCorrect,OpIds::PressureCorrectionCommit},
        "Structured HOW did not execute the multi-operation pressure fragments.");
    bool wrongOwnerRejected=false;
    SF::Run::OpRegistry wrongOwner;
    int wronglyExecuted=0;
    for (const auto& operation:SolvePlanner::requiredOperations(methodPlan))
        wrongOwner.bind(operation,"flow.conservative",[&] { ++wronglyExecuted; });
    try { SF::Run::PlanExecutor::execute(methodPlan,wrongOwner); }
    catch (const std::runtime_error& error) {
        wrongOwnerRejected=std::string(error.what()).find("frozen provider")!=std::string::npos;
    }
    require(wrongOwnerRejected && wronglyExecuted==0,
            "Runtime replaced the compiled provider or executed before ownership validation.");
    const auto frozenOwners=[&](const auto& self,const SolvePlanNode& node)->void {
        if (!node.operation.empty())
            require(node.provider==(node.operation==OpIds::FluxCorrect
                ? "flow.rhie-chow" : "flow.pressure-operators"),
                "Pressure method did not freeze leaf ownership, including lifecycle/publication.");
        for (const auto& child:node.children) self(self,child);
    };
    frozenOwners(frozenOwners,methodPlan.root);
    bool pressureRejected=false;
    try {
        ExecutionProgram missing{pressureProgram.root.children.at(1).step};
        (void)compileExecutionProgram(pressureSystem,missing,
            {{"pSimple","PressureCorrection"}},methods);
    } catch (const std::runtime_error&) { pressureRejected=true; }
    require(pressureRejected,"PressureCorrection accepted missing face-response producer.");
    pressureRejected=false;
    try {
        auto missingInput=pressureBindings;
        missingInput.front().inputs={"momentum"};
        (void)compileExecutionProgram(pressureSystem,pressureProgram,missingInput,methods);
    } catch (const std::runtime_error&) { pressureRejected=true; }
    require(pressureRejected,"PressureMomentum accepted duplicate source mathematics.");
    pressureRejected=false;
    try {
        auto wrongTarget=pressureProgram;
        wrongTarget.root.children.at(2).step.target.kind=TargetKind::Working;
        (void)compileExecutionProgram(pressureSystem,wrongTarget,pressureBindings,methods);
    } catch (const std::runtime_error&) { pressureRejected=true; }
    require(pressureRejected,"Pressure correction kernel accepted a target storage kind it cannot implement.");

    // Repeated calls share one equation definition while retaining independently
    // ordered storage and numerical bindings. The callbacks exercise only IR.
    ExecutionProgram repeated;
    ExecutionScope scope;
    scope.kind=ExecutionKind::Loop;
    scope.id="local";scope.order=10;scope.repetitions=2;
    ExecutionScope working;
    working.kind=ExecutionKind::EquationCall;
    working.order=20;
    working.step={"UpdateC",{"C",TargetKind::Working},"predict"};
    ExecutionScope physical;
    physical.kind=ExecutionKind::EquationCall;
    physical.order=60;
    physical.step={"UpdateC",{"C"},"publish"};
    scope.children={physical,working};
    repeated.root.children={scope};
    const auto mathBefore=formulaText(system.registry.at("UpdateC"));
    const auto repeatedCompiled=compileExecutionProgram(system,repeated,
        {{"*","MissingMethod"},{"UpdateC","DirectEvaluation"},
         {"UpdateC","LinearEquation",{},"publish"}},methods);
    require(repeatedCompiled.steps.size()==2
        && repeatedCompiled.steps[0].source.equation==repeatedCompiled.steps[1].source.equation
        && repeatedCompiled.steps[0].equationMethod=="DirectEvaluation"
        && repeatedCompiled.steps[1].equationMethod=="LinearEquation"
        && repeatedCompiled.steps[0].target.resources.front().storage=="state.C.Working"
        && repeatedCompiled.steps[1].target.resources.front().storage=="testStorage",
        "Repeated equation occurrences lost local order, override or typed storage.");
    require(repeated.root.children.front().children.front().step.occurrence=="publish"
        && system.registry.entries().size()==3 && formulaText(system.registry.at("UpdateC"))==mathBefore,
        "Compiler changed frozen source HOW/WHAT while sorting or binding storage.");
    const auto reorderedBindings=compileExecutionProgram(system,repeated,
        {{"UpdateC","LinearEquation",{},"publish"},{"UpdateC","DirectEvaluation"},
         {"*","MissingMethod"},{"*","AnotherMissingDefault"}},methods);
    require(reorderedBindings.steps[0].equationMethod==repeatedCompiled.steps[0].equationMethod
        && reorderedBindings.steps[1].equationMethod==repeatedCompiled.steps[1].equationMethod,
        "Highest-priority selection depends on binding registration order.");
    require(targetFromSyntax("U*").symbol=="U" && targetFromSyntax("U*").kind==TargetKind::Working
        && targetFromSyntax("p'").symbol=="p" && targetFromSyntax("p'").kind==TargetKind::Correction,
        "Source target notation was not converted to semantics before runtime.");
    CompiledSolvePlan repeatedPlan;
    repeatedPlan.root=compileMethodProgram(repeatedCompiled);
    std::vector<TargetKind> observedTargets;
    SF::Run::OpRegistry repeatedOperations;
    for (const auto& occurrence:repeatedCompiled.steps)
        repeatedOperations.bind(occurrence.backendOperation,
            [&](const SF::Run::ExecutionContext& context) {
                require(context.target!=nullptr,"Executor dropped a compiled target.");
                observedTargets.push_back(context.target->kind);
            });
    SF::Run::PlanExecutor::execute(repeatedPlan,repeatedOperations);
    require(observedTargets==std::vector<TargetKind>{TargetKind::Working,TargetKind::Physical,
            TargetKind::Working,TargetKind::Physical},
        "Generic loop did not preserve per-occurrence typed target execution.");
    bool rejected=false;
    try {
        (void)compileExecutionProgram(system,
            ExecutionProgram{EquationCall{"MyScalar",{""}}},binding,methods);
    } catch (const std::runtime_error&) { rejected=true; }
    require(rejected,"Compiler inferred a missing primary Output.");
    rejected=false;
    try {
        (void)compileExecutionProgram(system,program,
            {{"MyScalar","MissingMethod"}},methods);
    } catch (const std::runtime_error&) { rejected=true; }
    require(rejected,"Unknown numerical method silently fell back.");
    rejected=false;
    try {
        (void)compileExecutionProgram(system,program,
            {{"MyScalar","ConservativeResidual"},
             {"Unused","DirectEvaluation"}},methods);
    } catch (const std::runtime_error&) { rejected=true; }
    require(rejected,"Unused method selection was silently ignored.");
    ExecutionProgram structured;
    ExecutionScope repeat;
    repeat.kind=ExecutionKind::Loop;
    repeat.repetitions=2;
    ExecutionScope producer;
    producer.kind=ExecutionKind::EquationCall;
    producer.step={"UpdateC",{"C"}};
    repeat.children.push_back(producer);
    ExecutionScope consumer;
    consumer.kind=ExecutionKind::EquationCall;
    consumer.step={"DiffuseC",{"C"}};
    structured.root.children={repeat,consumer};
    const auto nested=compileExecutionProgram(system,structured,
        {{"UpdateC","DirectEvaluation"},{"DiffuseC","LinearEquation"}},methods);
    require(nested.root.children.size()==2
            && nested.root.children.front().kind==ExecutionKind::Loop
            && nested.steps.size()==2
,
            "Structured HOW did not preserve nesting and workspace dependency.");
    CompiledSolvePlan nestedPlan;
    nestedPlan.root=compileMethodProgram(nested);
    std::vector<std::string> order;
    SF::Run::OpRegistry operations;
    operations.bind(nested.steps.front().backendOperation,
        [&] { order.push_back("UpdateC"); });
    operations.bind(nested.steps.back().backendOperation,
        [&] { order.push_back("DiffuseC"); });
    SF::Run::PlanExecutor::execute(nestedPlan,operations);
    require(order==std::vector<std::string>{"UpdateC","UpdateC","DiffuseC"},
            "Structured HOW Repeat/Sequence did not execute in declared order.");
    ExecutionProgram until;
    ExecutionScope loop;
    loop.kind=ExecutionKind::Loop;
    loop.repetitions=3;
    loop.terminationSignal="converged";
    loop.children.push_back(producer);
    until.root.children.push_back(loop);
    const auto untilCompiled=compileExecutionProgram(system,until,
        {{"UpdateC","DirectEvaluation"}},methods);
    CompiledSolvePlan untilPlan;
    untilPlan.root=compileMethodProgram(untilCompiled);
    int iterations=0;
    SF::Run::OpRegistry untilOperations;
    untilOperations.bind(untilCompiled.steps.front().backendOperation,
        [&](const SF::Run::ExecutionContext& context) {
            ++iterations;
            context.signals->publish("converged",true);
        });
    SF::Run::PlanExecutor::execute(untilPlan,untilOperations);
    require(iterations==1,"HOW LoopUntil did not honor compiled convergence signal.");
    rejected=false;
    try {
        (void)compileExecutionProgram(system,structured,
            {{"UpdateC","DirectEvaluation"},{"UpdateC","DirectEvaluation"},
             {"DiffuseC","LinearEquation"}},methods);
    } catch (const std::runtime_error&) { rejected=true; }
    require(rejected,"Duplicate method binding was accepted.");
    for (const auto& id:std::vector<std::string>{"relaxIterate","restoreFlux","checkConvergence"})
        pressureSystem.registry.add({id,Expr::symbol(id),Expr::constantValue(0),{}});
    pressureSystem.state.add(stored("C"));
    pressureSystem.registry.add(system.registry.at("UpdateC"));
    for (const auto preset:{SF::FDM::PressureCouplingPreset::PISO,
            SF::FDM::PressureCouplingPreset::SIMPLE,SF::FDM::PressureCouplingPreset::PIMPLE}) {
        CouplingPresetRequest request;
        request.presetKind=preset;request.preset=SF::FDM::toString(preset);
        request.outerCorrectors=preset==SF::FDM::PressureCouplingPreset::PISO ? 1 : 3;
        request.pressureCorrectors=preset==SF::FDM::PressureCouplingPreset::SIMPLE ? 1 : 2;
        ExecutionProgram how;
        ExecutionScope defaultMomentum;
        defaultMomentum.kind=ExecutionKind::EquationCall;
        defaultMomentum.order=20;defaultMomentum.step={"momentum",{"U"}};
        ExecutionScope other;
        other.kind=ExecutionKind::EquationCall;other.order=70;other.step={"UpdateC",{"C"}};
        how.root.children={defaultMomentum,other};
        const auto registrySize=pressureSystem.registry.entries().size();
        const auto which=pressureNumerics(request);
        applyPressureExecution(how,request);
        require(pressureSystem.registry.entries().size()==registrySize
            && which.front().occurrence=="predictor",
            "Pressure HOW changed WHAT or occurrence numerical selection.");
        auto selected=which;
        selected.push_back({"UpdateC","DirectEvaluation"});
        const auto lowered=compileExecutionProgram(pressureSystem,how,selected,methods,&euler);
        CompiledSolvePlan plan;plan.root=compileMethodProgram(lowered);
        std::vector<std::string> actual;
        SF::Run::OpRegistry fake;
        int convergences=0;
        for (const auto& operation:SolvePlanner::requiredOperations(plan))
            fake.bind(operation,operation==OpIds::FluxCorrect ? "flow.rhie-chow"
                : operation.compare(0,8,"formula.")==0 ? "" : "flow.pressure-operators",[&,operation](const SF::Run::ExecutionContext& context) {
                actual.push_back(operation);
                if (operation==OpIds::PressureConvergenceEvaluate)
                    context.signals->publish(kPressureOuterConvergedSignal,++convergences==2);
            });
        SF::Run::PlanExecutor::execute(plan,fake);
        std::vector<std::string> expected{OpIds::PressurePrepare};
        const bool fixed=preset!=SF::FDM::PressureCouplingPreset::PISO;
        if (fixed) expected.push_back(OpIds::PressureStepBegin);
        for (int outer=0;outer<(fixed?2:1);++outer) {
            if (fixed) expected.push_back(OpIds::PressureIterationBegin);
            expected.insert(expected.end(),{OpIds::MomentumAssemble,OpIds::MomentumSolve});
            for (int correction=0;correction<request.pressureCorrectors;++correction)
                expected.insert(expected.end(),{OpIds::PressureBoundaryPrepare,OpIds::PressureAssemble,
                    OpIds::PressureSolve,OpIds::VelocityCorrect,OpIds::PressureUpdatePrepare,
                    OpIds::FluxCorrect,OpIds::PressureCorrectionCommit});
            if (fixed) expected.insert(expected.end(),{OpIds::PressureRelaxationApply,
                OpIds::PressureFluxConsistencyRestore,OpIds::PressureConvergenceEvaluate,
                OpIds::PressureIterationEnd});
        }
        expected.insert(expected.end(),{"formula.direct.evaluate.UpdateC",OpIds::PressureStepCommit,OpIds::TimeCommit});
        require(actual==expected,"Provider lifecycle changed pressure operation order/count/termination/commit.");
        bool temporalRejected=false;
        try {
            (void)compileExecutionProgram(pressureSystem,how,selected,methods,&rk4);
        } catch (const std::runtime_error& error) {
            temporalRejected=std::string(error.what()).find("Unsupported: provider")!=std::string::npos;
        }
        require(temporalRejected,"Unsupported pressure target/temporal capability silently fell back.");
    }
    ExecutionScope sourceCommit;
    sourceCommit.kind=ExecutionKind::Commit;
    sourceCommit.id="explicitCommit";
    auto unimplementedCommit=assignment;
    unimplementedCommit.root.children.push_back(sourceCommit);
    rejected=false;
    try {
        (void)compileExecutionProgram(system,unimplementedCommit,{{"UpdateC","DirectEvaluation"}},methods);
    } catch (const std::runtime_error&) { rejected=true; }
    require(rejected,"Source Commit without provider implementation was silently ignored.");
    ExecutionProgram nestedCommit;
    ExecutionScope nestedBody;
    nestedBody.children={program.root.children.front(),sourceCommit};
    nestedCommit.root.children={nestedBody};
    rejected=false;
    try {
        (void)compileExecutionProgram(system,nestedCommit,binding,methods,&euler,&temporalMethods.at(euler.id()));
    } catch (const std::runtime_error&) { rejected=true; }
    require(rejected,"Temporal fusion dropped a nested source Commit.");
    for (auto id:{SF::FDM::TimeRecipeId::ForwardEuler,
                 SF::FDM::TimeRecipeId::SSPRK3,
                 SF::FDM::TimeRecipeId::ClassicalRK4}) {
        const auto temporalMethods=builtinTemporalMethods();
        const auto& temporalMethod=temporalMethods.at(id);
        const auto recipe=temporalMethod.compile(
            SF::FDM::builtInTimeRecipe(id));
        auto temporalProgram=program;
        const auto temporal=compileExecutionProgram(system,temporalProgram,
            binding,methods,&recipe,&temporalMethod);
        require(temporal.hasTemporalRoot
                && temporal.steps.front().source.equation=="MyScalar"
                && temporal.temporalRoot.children.at(3).repetitions==recipe.stageCount(),
                "TemporalMethod did not wrap the same WHAT/HOW with its own stages.");
    }
    return 0;
}
