#include "solver/system/SF_methodObjects.h"
#include "solver/run/SF_planExecutor.h"
#include "core/system/SF_operationIds.h"
#include <stdexcept>
#include <vector>
using namespace SF::System;
namespace {
void require(bool condition,const char* message) { if (!condition) throw std::runtime_error(message); }
class Independent final : public IProvider {
public:
    explicit Independent(bool rootDecoration=false): rootDecoration_(rootDecoration) {}
    std::string_view id() const override { return "Independent"; }
    std::string_view runtimeProvider() const override { return "test.independent"; }
    void lowerLifecycle(const ExecutionScope&,bool root,const std::vector<CompiledEquationCall>&,
                        SolvePlanNode& node) const override {
        if (root && rootDecoration_) {
            SolvePlanNode extra;extra.kind=PlanNodeKind::Update;extra.operation="test.extra";
            extra.provider=runtimeProvider();node.children.push_back(extra);
        }
    }
    CompiledEquationCall compile(const ExecutableEquationSystem&,const EquationCall& call,
                                const NumericalBinding&) const override {
        CompiledEquationCall result;
        result.source=call;result.target.symbol=call.target.symbol;result.target.kind=call.target.kind;
        result.equationMethod=std::string(id());result.backendOperation="test.once";
        result.calls={{call.equation,call.target.symbol}};
        return result;
    }
private:
    bool rootDecoration_;
};
}
int main() {
    ExecutableEquationSystem system;
    for (const auto* name:{"A","B"}) {
        StateSymbol symbol;symbol.id=name;symbol.storageKey=name;
        system.state.add(symbol);
    }
    system.registry.add({"refresh",FormulaExpr::symbol("A"),FormulaExpr::constantValue(2),{}});
    system.registry.add({"advance",FormulaExpr::op("ddt",{FormulaExpr::symbol("B")}),FormulaExpr::constantValue(0),{}});
    ExecutionProgram source{{{"refresh",{"A"}},{"advance",{"B"}}}};
    ExecutionScope commit;commit.kind=ExecutionKind::Commit;source.root.children.push_back(commit);
    auto providers=builtinProviders();Independent independent;providers.add(independent);
    auto temporal=builtinTemporalMethods();const auto& method=temporal.at(SF::FDM::TimeRecipeId::ClassicalRK4);
    const auto recipe=method.compile(SF::FDM::builtInTimeRecipe(method.id()));
    const auto compile=[&](const ExecutionProgram& program) {
        return compileExecutionProgram(system,program,{{"refresh","Independent"},{"advance","ConservativeResidual"}},
            providers,&recipe,&method);
    };
    const auto compiled=compile(source);
    CompiledSolvePlan plan;plan.root=compiled.temporalRoot;
    SF::Run::OpRegistry callbacks;std::vector<std::string> order;
    for (const auto* op:{OpIds::FlowStepPrepare,OpIds::FlowDtCompute,OpIds::FlowStepBegin,
        OpIds::ExplicitStageExecute,OpIds::FlowStepCommit,OpIds::TimeCommit})
        callbacks.bind(op,"flow.conservative",[&,op] { order.push_back(op); });
    callbacks.bind("test.once","test.independent",[&] { order.push_back("test.once"); });
    SF::Run::PlanExecutor::execute(plan,callbacks);
    require(order==std::vector<std::string>{OpIds::FlowStepPrepare,OpIds::FlowDtCompute,"test.once",
        OpIds::FlowStepBegin,OpIds::ExplicitStageExecute,OpIds::ExplicitStageExecute,
        OpIds::ExplicitStageExecute,OpIds::ExplicitStageExecute,OpIds::FlowStepCommit,OpIds::TimeCommit},
        "Independent physical-step prefix was reordered, dropped or staged.");
    auto suffix=source;
    ExecutionScope group;group.kind=ExecutionKind::Sequence;group.id="post-predictor";
    auto after=source.root.children.front();after.step.occurrence="afterPredictor";
    group.children.push_back(after);
    after.step.occurrence="afterProjection";group.children.push_back(after);
    suffix.root.children.insert(suffix.root.children.end()-1,group);
    order.clear();plan.root=compile(suffix).temporalRoot;
    SF::Run::PlanExecutor::execute(plan,callbacks);
    require(order==std::vector<std::string>{OpIds::FlowStepPrepare,OpIds::FlowDtCompute,"test.once",
        OpIds::FlowStepBegin,OpIds::ExplicitStageExecute,OpIds::ExplicitStageExecute,
        OpIds::ExplicitStageExecute,OpIds::ExplicitStageExecute,"test.once","test.once",
        OpIds::FlowStepCommit,OpIds::TimeCommit},
        "Authored post-predictor group was staged, reordered or moved past physical commit.");
    auto interleaved=source;
    std::swap(interleaved.root.children[0],interleaved.root.children[1]);
    auto last=interleaved.root.children[0];last.step.occurrence="another";
    interleaved.root.children.insert(interleaved.root.children.end()-1,last);
    bool rejected=false;
    try { (void)compile(interleaved); } catch (const std::runtime_error& error) {
        rejected=std::string(error.what()).find("Unsupported mixed temporal topology")!=std::string::npos;
    }
    require(rejected,"Interleaved independent/temporal calls were silently reordered.");
    auto noCommit=source;noCommit.root.children.pop_back();rejected=false;
    try { (void)compile(noCommit); } catch (const std::runtime_error&) { rejected=true; }
    require(rejected,"Mixed physical-step composition accepted missing Commit.");
    auto decoratedProviders=builtinProviders();Independent decorated(true);decoratedProviders.add(decorated);
    rejected=false;
    try {
        (void)compileExecutionProgram(system,source,{{"refresh","Independent"},{"advance","ConservativeResidual"}},
            decoratedProviders,&recipe,&method);
    } catch (const std::runtime_error& error) {
        rejected=std::string(error.what()).find("root lifecycle decorations")!=std::string::npos;
    }
    require(rejected,"Mixed lowering silently discarded a provider root lifecycle operation.");
}
