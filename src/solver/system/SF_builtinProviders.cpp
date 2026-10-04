/// Numerical providers own domain matching, storage realization and local lifecycle.
#include "SF_methodObjects.h"
#include "SF_immersedMethods.h"
#include "SF_eulerianCoupling.h"
#include "SF_pressureCoupling.h"
#include "core/system/SF_operationIds.h"
#include <algorithm>
#include <stdexcept>
#include <utility>
namespace SF::System {
namespace {

SolvePlanNode methodLeaf(PlanNodeKind kind,std::string id,
                         std::string name,std::string operation,
                         std::string_view provider = {}) {
    SolvePlanNode leaf;
    leaf.kind=kind;
    leaf.id=std::move(id);
    leaf.name=std::move(name);
    leaf.operation=std::move(operation);
    leaf.provider=provider;
    return leaf;
}

SolvePlanNode explicitFragment(const CompiledEquationCall& body,
                               const CompiledTimeRecipe& recipe,
                               const std::vector<SolvePlanNode>& prefix) {
    if (!body.temporalResidual || body.calls.empty()
        || body.backendOperation.empty())
        throw std::runtime_error("TemporalMethod requires a compiled residual body and backend.");
    SolvePlanNode root;
    root.kind=PlanNodeKind::Sequence;
    root.id="Explicit.step";
    root.name="explicit time step";
    root.children.push_back(methodLeaf(PlanNodeKind::Update,
        "Explicit.prepare","prepare physical step",OpIds::FlowStepPrepare,body.backendProvider));
    root.children.push_back(methodLeaf(PlanNodeKind::Update,
        "Explicit.dt","compute stable time step",OpIds::FlowDtCompute,body.backendProvider));
    root.children.insert(root.children.end(),prefix.begin(),prefix.end());
    root.children.push_back(methodLeaf(PlanNodeKind::Update,
        "Explicit.begin","begin explicit step",OpIds::FlowStepBegin,body.backendProvider));
    auto stages=methodLeaf(PlanNodeKind::StageLoop,"Explicit.stages",
        std::string(FDM::toString(recipe.id()))+" fused explicit stages",{});
    stages.repetitions=recipe.stageCount();
    auto stage=methodLeaf(PlanNodeKind::Update,"Explicit.stage.execute",
        "execute explicit stage",body.backendOperation,body.backendProvider);
    stage.equationCalls=body.calls;
    stages.children.push_back(std::move(stage));
    root.children.push_back(std::move(stages));
    return root;
}

class ForwardEulerMethod final : public ITemporalMethod {
public:
    FDM::TimeRecipeId id() const override { return FDM::TimeRecipeId::ForwardEuler; }
    CompiledTimeRecipe compile(FDM::TimeRecipe selected) const override {
        return CompiledTimeRecipe::fromMethod(selected,{{{0.0,0.0,1.0}}},1,
            ExplicitStageBackend::ForwardEuler);
    }
    SolvePlanNode compileFragment(const CompiledTimeRecipe& compiled,
                                  const CompiledEquationCall& body,
                                  const std::vector<SolvePlanNode>& prefix) const override {
        return explicitFragment(body,compiled,prefix);
    }
};

class SSPRK3Method final : public ITemporalMethod {
public:
    FDM::TimeRecipeId id() const override { return FDM::TimeRecipeId::SSPRK3; }
    CompiledTimeRecipe compile(FDM::TimeRecipe selected) const override {
        return CompiledTimeRecipe::fromMethod(selected,
            {{{0.0,0.0,1.0},{1.0,0.75,0.25},{0.5,1.0/3.0,2.0/3.0}}},3,
            ExplicitStageBackend::SSPRK3);
    }
    SolvePlanNode compileFragment(const CompiledTimeRecipe& compiled,
                                  const CompiledEquationCall& body,
                                  const std::vector<SolvePlanNode>& prefix) const override {
        return explicitFragment(body,compiled,prefix);
    }
};

class ClassicalRK4Method final : public ITemporalMethod {
public:
    FDM::TimeRecipeId id() const override { return FDM::TimeRecipeId::ClassicalRK4; }
    CompiledTimeRecipe compile(FDM::TimeRecipe selected) const override {
        return CompiledTimeRecipe::fromMethod(selected,
            {{{0.0,0.0,0.5},{0.5,0.0,0.5},{0.5,0.0,1.0},{1.0,0.0,0.0}}},
            4,ExplicitStageBackend::ClassicalRK4,
            {{1.0,2.0,2.0,1.0}},6.0);
    }
    SolvePlanNode compileFragment(const CompiledTimeRecipe& compiled,
                                  const CompiledEquationCall& body,
                                  const std::vector<SolvePlanNode>& prefix) const override {
        return explicitFragment(body,compiled,prefix);
    }
};

bool containsDdt(const FormulaExpr& expression, std::string_view target) {
    if (expression.kind==FormulaExpr::Kind::Operator
        && expression.name=="ddt" && expression.arguments.size()==1
        && expression.arguments.front().kind==FormulaExpr::Kind::Symbol
        && expression.arguments.front().name==target) return true;
    return std::any_of(expression.arguments.begin(),expression.arguments.end(),
        [&](const FormulaExpr& child) { return containsDdt(child,target); });
}

void collectSymbols(const FormulaExpr& expression,std::vector<std::string>& result) {
    if (expression.kind==FormulaExpr::Kind::Symbol
        && std::find(result.begin(),result.end(),expression.name)==result.end())
        result.push_back(expression.name);
    for (const auto& child:expression.arguments) collectSymbols(child,result);
}

void requireOutput(const ExecutableEquationSystem&,const EquationCall& step) {
    if (step.equation.empty() || step.target.symbol.empty())
        throw std::runtime_error("EquationCall requires an equation and a STATE target.");
}

CompiledTarget bindTarget(const ExecutableEquationSystem&,const EquationCall& call,
                          std::string workspace = {}) {
    CompiledTarget result;
    result.symbol=call.target.symbol;result.kind=call.target.kind;
    if (call.target.kind!=TargetKind::Physical)
        result.workspace=workspace.empty()
            ? "state."+call.target.symbol+"."+(call.target.kind==TargetKind::Working
                ? "Working" : call.target.kind==TargetKind::Correction ? "Correction" : "Workspace")
            : std::move(workspace);
    return result;
}

class ConservativeResidualMethod final : public IProvider {
public:
    bool usesSpatialRecipes() const override { return true; }
    std::string_view id() const override { return "ConservativeResidual"; }
    std::string_view runtimeProvider() const override { return "flow.conservative"; }
    void lowerLifecycle(const ExecutionScope& scope,bool,
            const std::vector<CompiledEquationCall>&,SolvePlanNode& node) const override {
        if (scope.kind!=ExecutionKind::Commit) return;
        node.children.push_back(methodLeaf(PlanNodeKind::Commit,"Explicit.commit",
            "commit physical state",OpIds::FlowStepCommit,runtimeProvider()));
        node.children.push_back(methodLeaf(PlanNodeKind::Commit,"Explicit.time.commit",
            "commit physical clock",OpIds::TimeCommit,runtimeProvider()));
    }
    CompiledEquationCall compile(const ExecutableEquationSystem& system,
                                const EquationCall& step,
                                const NumericalBinding&) const override {
        requireOutput(system,step);
        if (step.target.kind!=TargetKind::Physical)
            throw std::runtime_error("Unsupported: conservative residual provider requires a Physical target.");
        if (step.target.symbol.empty())
            throw std::runtime_error("ConservativeResidual requires one equation and one target.");
        CompiledEquationCall result;
        result.source=step;
        result.target=bindTarget(system,step);
        result.equationMethod=std::string(id());
        result.temporalResidual=true;
        result.oldTimeWorkspace="explicit.q0";
        result.publishesStageToPhysicalTarget=true;
        result.spatialTerms=true;
        result.residualWorkspace="residual";
        for (const char* operation:{OpIds::FlowStepPrepare,OpIds::FlowDtCompute,OpIds::FlowStepBegin,
                OpIds::ExplicitStageExecute,OpIds::FlowStepCommit,OpIds::TimeCommit})
            result.operations.push_back({operation,operation,OperationStage::Prepare,
                {OperationCapability::ConservativeExplicit},{OriginKind::Generated,"ConservativeResidual"}});
        result.backendOperation=OpIds::ExplicitStageExecute;
        result.writes={step.target.symbol};
        result.requirements={"boundary-ready state","conservative residual"};
        const auto& formula=system.registry.at(step.equation);
        const auto& symbol=step.target.symbol;
        if (!containsDdt(formula.lhs,symbol))
            throw std::runtime_error("ConservativeResidual equation has no ddt(target).");
        result.calls.push_back({formula.id,symbol});
        collectSymbols(formula.lhs,result.reads);
        collectSymbols(formula.rhs,result.reads);
        return result;
    }
};

// Capability matcher for the frozen backend, not an alternate WHAT registry.
Equation supportedRasMathematics(const std::string& second,
                                            const std::string& target) {
    using E=FormulaExpr;
    const auto s=[](const std::string& name) { return E::symbol(name); };
    const auto mul=[](E a,E b) { return E::multiply(std::move(a),std::move(b)); };
    const auto div=[](E a,E b) { return E::divide(std::move(a),std::move(b)); };
    if ((second!="epsilon" && second!="omega") || (target!="k" && target!=second))
        throw std::runtime_error("Unsupported single-fluid RAS mathematical contract.");
    const bool epsilon=second=="epsilon";
    const std::string model=epsilon ? "kEpsilon" : "kOmegaSST";
    const auto coefficient=[&](const std::string& name) { return s(model+"."+name); };
    const auto production=coefficient("limitedProduction");
    E source;
    E diffusivity;
    if (target=="k") {
        const auto sink=epsilon ? s("epsilon")
            : mul(coefficient("betaStar"),mul(s("k"),s("omega")));
        source=E::subtract(div(production,s("rho")),sink);
    } else if (epsilon) {
        source=E::subtract(div(mul(coefficient("C1"),mul(s("epsilon"),production)),
                                  coefficient("maxRhoK")),
            div(mul(coefficient("C2"),mul(s("epsilon"),s("epsilon"))),coefficient("maxK")));
    } else {
        source=E::subtract(div(mul(coefficient("gamma"),production),coefficient("maxMuT")),
            mul(coefficient("beta"),mul(s("omega"),s("omega"))));
    }
    diffusivity=epsilon ? div(s("mu_t"),coefficient(target=="k" ? "sigmaK" : "sigmaEpsilon"))
        : E::add(coefficient("laminarMu"),mul(coefficient(target=="k" ? "sigmaK" : "sigmaOmega"),s("mu_t")));
    auto diffusion=E::op("diffusion",{diffusivity,s(target)});
    if (!epsilon && target=="omega")
        diffusion=E::add(std::move(diffusion),coefficient("crossDiffusion"));
    return {target,E::op("ddt",{s(target)}),
        E::add(std::move(source),div(std::move(diffusion),s("rho"))),
        {OriginKind::Model,model},true};
}

class TurbulenceTransportMethod final : public IProvider {
public:
    std::string_view id() const override { return "TurbulenceTransport"; }
    std::string_view runtimeProvider() const override { return "flow.turbulence"; }
    CompiledEquationCall compile(const ExecutableEquationSystem& system,
                                const EquationCall& step,
                                const NumericalBinding& binding) const override {
        if (binding.inputs.size()!=2 || binding.inputs[0]!="k"
            || (binding.inputs[1]!="omega" && binding.inputs[1]!="epsilon"))
            throw std::runtime_error("TurbulenceTransport requires an explicit complete RAS pair.");
        const auto& second=binding.inputs[1];
        const auto model=second=="omega" ? "kOmegaSST" : "kEpsilon";
        if (system.state.contains(second=="omega" ? "epsilon" : "omega"))
            throw std::runtime_error("TurbulenceTransport has incompatible model/state pair.");
        if (step.target.kind!=TargetKind::Physical || step.target.symbol!=step.equation
            || (step.equation!="k" && step.equation!=second))
            throw std::runtime_error("TurbulenceTransport requires each equation's own Physical target.");
        for (const auto& name:binding.inputs) {
            const auto& state=system.state.at(name);
            if (state.role!=StateRole::Transported || state.storageBinding!=StorageBinding::ProviderDistributed
                || state.storageKey!=name || state.components!=1)
                throw std::runtime_error("TurbulenceTransport requires provider-owned transported STATE: "+name);
            if (canonicalFormula(system.registry.at(name))
                !=canonicalFormula(supportedRasMathematics(second,name)))
                throw std::runtime_error("Unsupported TurbulenceTransport mathematics for "+name
                    +"; the frozen kernel implements source/diffusion with no convection.");
        }
        const auto& closure=system.state.at("mu_t");
        if (closure.role!=StateRole::Derived || closure.storageBinding!=StorageBinding::ProviderDistributed
            || closure.storageKey!="mu_t")
            throw std::runtime_error("TurbulenceTransport requires the provider-owned mu_t closure.");
        CompiledEquationCall result;
        result.source=step;result.target=bindTarget(system,step);
        result.equationMethod=std::string(id());
        result.backendOperation=OpIds::TurbulenceAdvance;
        result.calls={{step.equation,step.target.symbol}};
        result.writes={step.target.symbol,"mu_t"};
        result.reads={"rho","U","k",second,"mu_t"};
        result.fusionKey=std::string(id())+"/"+model;
        result.fusionMembers={{"k","k"},{second,second}};
        result.temporalMethod="physical-step explicit in-place update";
        result.requirements={"single-fluid serial single patch","complete adjacent RAS pair",
            "boundary -> WriteOwned -> ReadHalo -> correction -> boundary -> WriteOwned -> ReadHalo"};
        result.operations.push_back({OpIds::TurbulenceAdvance,"advance frozen RAS pair once",
            OperationStage::Prepare,{OperationCapability::SingleFluidTurbulenceTransport},
            {OriginKind::Generated,std::string(id())}});
        return result;
    }
};

// Smagorinsky refresh is an algebraic closure, never a transported RAS pair.
class TurbulenceClosureMethod final : public IProvider {
public:
    std::string_view id() const override { return "TurbulenceClosure"; }
    std::string_view runtimeProvider() const override { return "flow.turbulence-closure"; }
    CompiledEquationCall compile(const ExecutableEquationSystem& system,
                                const EquationCall& step,const NumericalBinding&) const override {
        const Equation expected{"mu_t",FormulaExpr::symbol("mu_t"),
            FormulaExpr::symbol("Smagorinsky.eddyDynamicViscosity"),{},true};
        if (step.equation!="mu_t" || step.target.symbol!="mu_t" || step.target.kind!=TargetKind::Physical
            || canonicalFormula(system.registry.at(step.equation))!=canonicalFormula(expected)
            || system.state.at("mu_t").storageBinding!=StorageBinding::ProviderDistributed)
            throw std::runtime_error("Unsupported TurbulenceClosure equation/state contract.");
        CompiledEquationCall result;
        result.source=step;result.target=bindTarget(system,step);
        result.equationMethod=std::string(id());
        result.backendOperation=OpIds::TurbulenceClosureRefresh;
        result.calls={{"mu_t","mu_t"}};result.writes={"mu_t"};result.reads={"rho","U"};
        result.operations.push_back({OpIds::TurbulenceClosureRefresh,"refresh algebraic eddy viscosity",
            OperationStage::Prepare,{OperationCapability::SingleFluidTurbulenceClosure},{}});
        return result;
    }
};

class DirectEvaluationMethod final : public IProvider {
public:
    std::string_view id() const override { return "DirectEvaluation"; }
    CompiledEquationCall compile(const ExecutableEquationSystem& system,
                                const EquationCall& step,
                                const NumericalBinding&) const override {
        requireOutput(system,step);
        if (step.target.symbol.empty())
            throw std::runtime_error("DirectEvaluation requires one primary Output symbol.");
        const auto& formula=system.registry.at(step.equation);
        const auto& symbol=step.target.symbol;
        if (formula.lhs.kind!=FormulaExpr::Kind::Symbol
            || formula.lhs.name!=symbol)
            throw std::runtime_error("DirectEvaluation Equation '"+formula.id
                +"' does not assign the declared primary Output '"+symbol+"'.");
        CompiledEquationCall result;
        result.source=step;
        result.target=bindTarget(system,step);
        result.equationMethod=std::string(id());
        result.backendOperation="formula.direct.evaluate."+step.equation;
        result.calls={{formula.id,symbol}};
        result.writes={symbol};
        if (!result.target.workspace.empty())
            result.workspaceProvides={result.target.workspace};
        collectSymbols(formula.rhs,result.reads);
        return result;
    }
};

class LinearEquationMethod final : public IProvider {
public:
    std::string_view id() const override { return "LinearEquation"; }
    CompiledEquationCall compile(const ExecutableEquationSystem& system,
                                const EquationCall& step,
                                const NumericalBinding&) const override {
        requireOutput(system,step);
        if (step.target.symbol.empty())
            throw std::runtime_error("LinearEquation requires one primary Output symbol.");
        const auto& formula=system.registry.at(step.equation);
        const auto& symbol=step.target.symbol;
        if (!dependsOn(formula.lhs,symbol)
            && !dependsOn(formula.rhs,symbol))
            throw std::runtime_error("LinearEquation Equation '"+formula.id
                +"' does not contain its declared Output '"+symbol+"'.");
        CompiledEquationCall result;
        result.source=step;
        result.target=bindTarget(system,step);
        result.equationMethod=std::string(id());
        result.backendOperation="formula.linear.solve."+step.equation;
        result.calls={{formula.id,symbol}};
        result.writes={symbol};
        if (!result.target.workspace.empty())
            result.workspaceProvides={result.target.workspace};
        result.requirements={"linear assembly provider","GlobalDofSystem",
                             "linear solver backend"};
        collectSymbols(formula.lhs,result.reads);
        collectSymbols(formula.rhs,result.reads);
        return result;
    }
};

/// Pressure methods compile only the local numerical topology. The existing
/// PressureOperators callbacks retain gauge, boundary and canonical-face math.
std::vector<ExecutableOperation> pressureOperations() {
    std::vector<ExecutableOperation> operations;
        const Provenance source{OriginKind::Generated,"pressureConstraint"};
        const auto declare = [&](OperationStage stage, const char* id,
                                 const char* name,
                                 std::vector<OperationCapability> needs) {
            operations.push_back(
                {id,name,stage,std::move(needs),source});
        };
        declare(OperationStage::Prepare,OpIds::PressurePrepare,
                "prepare pressure schedule",{OperationCapability::PressureSchedule});
        declare(OperationStage::FixedTimeStepBegin,OpIds::PressureStepBegin,
                "freeze physical-time base state",
                {OperationCapability::FixedTimeIteration});
        declare(OperationStage::IterationBegin,OpIds::PressureIterationBegin,
                "begin fixed-time outer iterate",
                {OperationCapability::FixedTimeIteration});
        declare(OperationStage::MomentumAssemble,OpIds::MomentumAssemble,
                "assemble momentum predictor",{OperationCapability::MomentumPredictor});
        declare(OperationStage::MomentumSolve,OpIds::MomentumSolve,
                "advance momentum predictor",{OperationCapability::MomentumPredictor});
        declare(OperationStage::PressureBoundaryPrepare,
                OpIds::PressureBoundaryPrepare,"prepare pressure boundary state",
                {OperationCapability::PressureBoundary});
        declare(OperationStage::PressureAssemble,OpIds::PressureAssemble,
                "assemble pressure correction",{OperationCapability::PressureCorrection,
                                                  OperationCapability::PressureLinearSolve});
        declare(OperationStage::PressureSolve,OpIds::PressureSolve,
                "solve pressure correction",{OperationCapability::PressureCorrection,
                                               OperationCapability::PressureLinearSolve});
        declare(OperationStage::PressureUpdatePrepare,
                OpIds::PressureUpdatePrepare,"prepare pressure update",
                {OperationCapability::PressureCorrection});
        declare(OperationStage::VelocityCorrect,OpIds::VelocityCorrect,
                "velocity correction",{OperationCapability::VelocityCorrection});
        declare(OperationStage::FluxCorrect,OpIds::FluxCorrect,
                "refresh derived face-flux state",{OperationCapability::FluxCorrection});
        declare(OperationStage::CorrectionCommit,
                OpIds::PressureCorrectionCommit,
                "commit pressure-corrected state",{OperationCapability::PressureCorrection});
        declare(OperationStage::RelaxationApply,OpIds::PressureRelaxationApply,
                "relax fixed-point solution",
                {OperationCapability::FixedTimeIteration});
        declare(OperationStage::FluxConsistencyRestore,
                OpIds::PressureFluxConsistencyRestore,
                "restore face flux for the relaxed iterate",
                {OperationCapability::FixedTimeIteration,
                 OperationCapability::FluxCorrection});
        declare(OperationStage::ConvergenceEvaluate,
                OpIds::PressureConvergenceEvaluate,
                "evaluate fixed-point residual",
                {OperationCapability::FixedTimeIteration});
        declare(OperationStage::IterationEnd,OpIds::PressureIterationEnd,
                "end fixed-time outer iterate",
                {OperationCapability::FixedTimeIteration});
        declare(OperationStage::StepCommit,OpIds::PressureStepCommit,
                "commit corrected state",{OperationCapability::PressureSchedule});
    operations.push_back({OpIds::TimeCommit,"commit physical clock",OperationStage::StepCommit,{OperationCapability::PressureSchedule},source});
    return operations;
}

class PressureMethod final : public IProvider {
public:
    bool usesSpatialRecipes() const override { return true; }
    enum class OutputContract { Assignment, Relation };
    PressureMethod(std::string_view method,std::string_view output,
                   OutputContract contract,std::vector<SolvePlanNode> leaves,
                   std::vector<std::string> needs,
                   std::vector<std::string> provides,
                   std::vector<std::string> requirements,
                   std::size_t sourceInputs=0)
        : id_(method),output_(output),contract_(contract),leaves_(std::move(leaves)),
          needs_(std::move(needs)),provides_(std::move(provides)),
          requirements_(std::move(requirements)),sourceInputs_(sourceInputs) {}
    std::string_view id() const override { return id_; }
    std::string_view runtimeProvider() const override {
        return id_=="FluxCorrection" ? "flow.rhie-chow" : "flow.pressure-operators";
    }
    std::string_view operationProvider(const OpId& operation) const override {
        // Flux publication uses the pressure backend's existing state authority.
        return operation==OpIds::FluxCorrect ? "flow.rhie-chow" : "flow.pressure-operators";
    }
    void lowerLifecycle(const ExecutionScope& scope,bool root,
            const std::vector<CompiledEquationCall>& calls,SolvePlanNode& node) const override {
        if (id_!="PressureMomentum") return;
        const bool fixed=std::any_of(calls.begin(),calls.end(),[](const auto& call) {
            return call.equationMethod=="FixedTimeRelaxation";
        });
        const auto leaf=[&](const char* operation,PlanNodeKind kind=PlanNodeKind::Update) {
            return methodLeaf(kind,node.id+"."+operation,operation,operation,operationProvider(operation));
        };
        if (root) {
            std::vector<SolvePlanNode> prefix{leaf(OpIds::PressurePrepare)};
            if (fixed) prefix.push_back(leaf(OpIds::PressureStepBegin));
            node.children.insert(node.children.begin(),prefix.begin(),prefix.end());
        }
        const auto ownsPredictor=[&](const auto& self,const ExecutionScope& source)->bool {
            if (source.kind==ExecutionKind::EquationCall)
                return std::any_of(calls.begin(),calls.end(),[&](const auto& call) {
                    return call.source.occurrence==source.step.occurrence && call.equationMethod==id_;
                });
            return std::any_of(source.children.begin(),source.children.end(),[&](const auto& child) { return self(self,child); });
        };
        if (fixed && scope.kind==ExecutionKind::Loop && !scope.terminationSignal.empty()
            && ownsPredictor(ownsPredictor,scope)) {
            node.children.insert(node.children.begin(),leaf(OpIds::PressureIterationBegin));
            node.children.push_back(leaf(OpIds::PressureIterationEnd));
        }
        if (scope.kind==ExecutionKind::Commit) {
            node.children.push_back(leaf(OpIds::PressureStepCommit,PlanNodeKind::Commit));
            node.children.push_back(leaf(OpIds::TimeCommit,PlanNodeKind::Commit));
        }
    }
    CompiledEquationCall compile(const ExecutableEquationSystem& system,
                                const EquationCall& step,
                                const NumericalBinding& binding) const override {
        requireOutput(system,step);
        const bool supportedTarget=id_=="PressureMomentum"
            ? (step.target.kind==TargetKind::Physical || step.target.kind==TargetKind::Working)
            : step.target.kind==(id_=="PressureCorrection" ? TargetKind::Correction
                : TargetKind::Physical);
        if (!supportedTarget)
            throw std::runtime_error("Unsupported: provider "+std::string(id_)+" cannot realize semantic target "+targetText(step.target));
        if (step.target.symbol.empty() || step.target.symbol!=
            (id_=="PressureCorrection" && step.target.kind==TargetKind::Correction ? "p" : output_))
            throw std::runtime_error("EquationMethod '"+std::string(id_)
                +"' requires explicit Output '"+std::string(output_)+"'.");
        const auto& formula=system.registry.at(step.equation);
        if (contract_==OutputContract::Relation) {
            if (!dependsOn(formula.lhs,output_)
                && !dependsOn(formula.rhs,output_))
                throw std::runtime_error("Pressure method Equation does not contain its Output.");
        } else if (formula.lhs.kind!=FormulaExpr::Kind::Symbol
                   || formula.lhs.name!=output_) {
            throw std::runtime_error("Pressure method Equation does not assign its Output.");
        }
        if (binding.inputs.size()!=sourceInputs_)
            throw std::runtime_error("EquationMethod '"+std::string(id_)
                +"' has missing or unexpected source-math inputs.");
        for (const auto& input:binding.inputs) {
            if (std::none_of(system.registry.entries().begin(),
                    system.registry.entries().end(),
                    [&](const Equation& candidate) { return candidate.id==input; }))
                throw std::runtime_error("EquationMethod source mathematics '"
                    +input+"' is unavailable.");
        }
        CompiledEquationCall result;
        result.source=step;
        result.target=bindTarget(system,step);
        result.equationMethod=std::string(id_);
        result.target=bindTarget(system,step,provides_.empty() ? "" : provides_.back());
        result.target.viewOwner=StateViewOwner::NumericalProvider;
        result.temporalCapabilities={FDM::TimeRecipeId::ForwardEuler};
        if (id_=="PressureMomentum") {
            result.spatialTerms=true;
            result.primitiveSourceRequired=true;
            result.residualWorkspace="pressureMomentumWorkspace";
            result.operations=pressureOperations();
        }
        result.calls={{formula.id,std::string(output_)}};
        result.writes={std::string(output_)};
        result.requirements=requirements_;
        result.workspaceRequires=needs_;
        result.workspaceProvides=provides_;
        result.sourceMathInputs=binding.inputs;
        result.reads=binding.inputs;
        collectSymbols(formula.lhs,result.reads);
        collectSymbols(formula.rhs,result.reads);
        result.fragment.kind=PlanNodeKind::Sequence;
        result.fragment.id=std::string(id_)+".fragment";
        result.fragment.name=std::string(id_)+" numerical fragment";
        result.fragment.children=leaves_;
        result.fragment.children.back().equationCalls=result.calls;
        return result;
    }
private:
    std::string_view id_;
    std::string_view output_;
    OutputContract contract_;
    std::vector<SolvePlanNode> leaves_;
    std::vector<std::string> needs_;
    std::vector<std::string> provides_;
    std::vector<std::string> requirements_;
    std::size_t sourceInputs_;
};

// The conservative predictor publishes the packed physical state in the
// frozen algorithm. Its fusion contract explicitly includes mass and energy;
// Working(p) and Correction(p) alias the Corrector's existing workspace.
class ConservativePressureMethod final : public IProvider {
public:
    bool usesSpatialRecipes() const override { return true; }
    ConservativePressureMethod(std::string_view method,std::string_view equation,
            std::string_view output,TargetKind kind,std::vector<const char*> operations,
            std::vector<std::string> needs={},std::vector<std::string> provides={})
        : method_(method),equation_(equation),output_(output),kind_(kind),
          operations_(std::move(operations)),needs_(std::move(needs)),provides_(std::move(provides)) {}
    std::string_view id() const override { return method_; }
    std::string_view runtimeProvider() const override { return "flow.conservative"; }
    void lowerLifecycle(const ExecutionScope& scope,bool root,
            const std::vector<CompiledEquationCall>& calls,SolvePlanNode& node) const override {
        if (method_!="ConservativePressureMomentum") return;
        const bool fixed=std::any_of(calls.begin(),calls.end(),[](const auto& call) {
            return call.equationMethod=="ConservativeFixedTimeRelaxation";
        });
        const auto leaf=[&](const char* operation,PlanNodeKind kind=PlanNodeKind::Update) {
            return methodLeaf(kind,node.id+"."+operation,operation,operation,runtimeProvider());
        };
        if (root) {
            std::vector<SolvePlanNode> prefix{leaf(OpIds::PressurePrepare)};
            if (fixed) prefix.push_back(leaf(OpIds::PressureStepBegin));
            node.children.insert(node.children.begin(),prefix.begin(),prefix.end());
        }
        const auto ownsPredictor=[&](const auto& self,const ExecutionScope& source)->bool {
            if (source.kind==ExecutionKind::EquationCall)
                return std::any_of(calls.begin(),calls.end(),[&](const auto& call) {
                    return call.source.occurrence==source.step.occurrence && call.equationMethod==method_;
                });
            return std::any_of(source.children.begin(),source.children.end(),[&](const auto& child) { return self(self,child); });
        };
        if (fixed && scope.kind==ExecutionKind::Loop && !scope.terminationSignal.empty()
            && ownsPredictor(ownsPredictor,scope)) {
            node.children.insert(node.children.begin(),leaf(OpIds::PressureIterationBegin));
            node.children.push_back(leaf(OpIds::PressureIterationEnd));
        }
        if (scope.kind==ExecutionKind::Commit) {
            node.children.push_back(leaf(OpIds::PressureStepCommit,PlanNodeKind::Commit));
            node.children.push_back(leaf(OpIds::TimeCommit,PlanNodeKind::Commit));
        }
    }
    CompiledEquationCall compile(const ExecutableEquationSystem& system,
            const EquationCall& call,const NumericalBinding& binding) const override {
        requireOutput(system,call);
        if (call.equation!=equation_ || call.target.symbol!=output_ || call.target.kind!=kind_)
            throw std::runtime_error("Unsupported: provider "+std::string(method_)+" requires "
                +std::string(equation_)+" -> "+std::string(output_)+" with its declared semantic target.");
        CompiledEquationCall result;
        result.source=call;result.equationMethod=method_;
        result.target=bindTarget(system,call,provides_.empty() ? "" : provides_.back());
        result.target.viewOwner=StateViewOwner::NumericalProvider;
        result.temporalCapabilities={FDM::TimeRecipeId::ForwardEuler};
        result.workspaceRequires=needs_;result.workspaceProvides=provides_;
        const auto& formula=system.registry.at(call.equation);
        if (method_=="ConservativePressureMomentum") {
            if (binding.inputs!=std::vector<std::string>{"continuity","energy"})
                throw std::runtime_error("ConservativePressureMomentum requires explicit continuity/energy fusion inputs.");
            result.spatialTerms=true;result.residualWorkspace="residual";
            result.operations=pressureOperations();
            result.sourceMathInputs=binding.inputs;
            for (const auto& math:std::vector<CompiledMathRef>{{"continuity","rho"},{"momentum","rhoU"},{"energy","rhoE"}}) {
                const auto& input=system.registry.at(math.equation);
                if (!containsDdt(input.lhs,math.target))
                    throw std::runtime_error("Conservative pressure predictor requires ddt("+math.target+").");
                result.calls.push_back(math);result.writes.push_back(math.target);
                collectSymbols(input.lhs,result.reads);collectSymbols(input.rhs,result.reads);
            }
            result.requirements={"five-component packed conservative state","EOS closure",
                "boundary-ready state","single-stage fused conservative residual"};
        } else {
            if (!binding.inputs.empty()) throw std::runtime_error("Unexpected conservative pressure source-math inputs.");
            const auto relations=conservativePressureRelations();
            const auto expected=std::find_if(relations.begin(),relations.end(),[&](const auto& value) { return value.id==call.equation; });
            if (expected==relations.end() || canonicalFormula(formula)!=canonicalFormula(*expected))
                throw std::runtime_error("Unsupported: selected conservative pressure provider cannot execute changed mathematics: "+call.equation);
            result.calls={{formula.id,std::string(output_)}};result.writes={std::string(output_)};
            collectSymbols(formula.lhs,result.reads);collectSymbols(formula.rhs,result.reads);
            result.requirements={"bound conservative pressure relation","EOS closure","pressure gauge/reference"};
        }
        result.fragment.kind=PlanNodeKind::Sequence;
        result.fragment.id=std::string(method_)+".fragment";
        result.fragment.name=std::string(method_)+" numerical fragment";
        for (const auto* operation:operations_)
            result.fragment.children.push_back(methodLeaf(PlanNodeKind::Update,
                std::string(method_)+"."+operation,operation,operation,runtimeProvider()));
        result.fragment.children.back().equationCalls=result.calls;
        return result;
    }
private:
    std::string_view method_,equation_,output_;
    TargetKind kind_;
    std::vector<const char*> operations_;
    std::vector<std::string> needs_,provides_;
};

class RelationProvider final : public IProvider {
public:
    RelationProvider(std::string id,std::string operation)
        : id_(std::move(id)),operation_(std::move(operation)) {}
    std::string_view id() const override { return id_; }
    std::string_view runtimeProvider() const override { return "flow.pressure-operators"; }
    CompiledEquationCall compile(const ExecutableEquationSystem& system,
            const EquationCall& call,const NumericalBinding&) const override {
        requireOutput(system,call);
        const auto& equation=system.registry.at(call.equation);
        CompiledEquationCall result;
        result.source=call; result.equationMethod=id_;
        result.target=bindTarget(system,call);
        result.target.viewOwner=StateViewOwner::NumericalProvider;
        result.temporalCapabilities={FDM::TimeRecipeId::ForwardEuler};
        result.backendOperation=operation_;
        result.calls.push_back({equation.id,call.target.symbol});
        result.writes.push_back(call.target.symbol);
        if (!result.target.workspace.empty()) result.workspaceProvides.push_back(result.target.workspace);
        collectSymbols(equation.lhs,result.reads);
        collectSymbols(equation.rhs,result.reads);
        return result;
    }
private:
    std::string id_,operation_;
};

} // namespace

TemporalMethodRegistry builtinTemporalMethods() {
    static const ForwardEulerMethod euler;
    static const SSPRK3Method ssp;
    static const ClassicalRK4Method rk4;
    TemporalMethodRegistry result;
    result.add(euler);
    result.add(ssp);
    result.add(rk4);
    return result;
}

ProviderRegistry builtinProviders() {
    static const ConservativeResidualMethod residual;
    static const DirectEvaluationMethod evaluation;
    static const LinearEquationMethod linear;
    static const PressureMethod momentum("PressureMomentum","U",PressureMethod::OutputContract::Relation,
        {methodLeaf(PlanNodeKind::Assemble,"PressureMomentum.assemble",
            "assemble momentum predictor",OpIds::MomentumAssemble),
         methodLeaf(PlanNodeKind::Solve,"PressureMomentum.solve",
            "solve momentum predictor",OpIds::MomentumSolve)},
        {},{"HbyA","rAU","faceResponse","correctedFlux","pressureMomentumWorkspace"},
        {"bound Equation occurrence mathematics",
         "velocity boundary","RhieChow face prediction","canonical face COPY"});
    static const PressureMethod correction("PressureCorrection","pPrime",
        PressureMethod::OutputContract::Relation,
        {methodLeaf(PlanNodeKind::Update,"PressureCorrection.boundary",
            "prepare pressure boundary state",OpIds::PressureBoundaryPrepare),
         methodLeaf(PlanNodeKind::Assemble,"PressureCorrection.assemble",
            "assemble pressure correction",OpIds::PressureAssemble),
         methodLeaf(PlanNodeKind::Solve,"PressureCorrection.solve",
            "solve pressure correction",OpIds::PressureSolve)},
        {"faceResponse","correctedFlux"},{"pressureCorrection"},
        {"pressure boundary","pressure gauge/reference","GlobalDofSystem",
         "linear solver backend"});
    static const PressureMethod pressureUpdate("PressureUpdate","p",
        PressureMethod::OutputContract::Assignment,
        {methodLeaf(PlanNodeKind::Correct,"PressureUpdate.prepare",
            "prepare pressure update",OpIds::PressureUpdatePrepare)},
        {"pressureCorrection"},{},
        {"pressure gauge/reference","pressure boundary"});
    static const PressureMethod velocityCorrection("VelocityCorrection","U",
        PressureMethod::OutputContract::Assignment,
        {methodLeaf(PlanNodeKind::Correct,"VelocityCorrection.correct",
            "correct velocity",OpIds::VelocityCorrect)},
        {"pressureCorrection","rAU"},{},
        {"velocity boundary closure"});
    static const PressureMethod fluxCorrection("FluxCorrection","phi",
        PressureMethod::OutputContract::Assignment,
        {methodLeaf(PlanNodeKind::Correct,"FluxCorrection.correct",
            "correct face flux",OpIds::FluxCorrect),
         methodLeaf(PlanNodeKind::Update,"FluxCorrection.publish",
            "publish corrected pressure state",OpIds::PressureCorrectionCommit)},
        {"pressureCorrection","faceResponse"},{"pressureFaceFlux"},
        {"RhieChow pressureFlux provider","canonical face COPY",
         "continuity defect"});
    static const RelationProvider relax("FixedTimeRelaxation",OpIds::PressureRelaxationApply);
    static const RelationProvider restore("FluxConsistency",OpIds::PressureFluxConsistencyRestore);
    static const RelationProvider convergence("ResidualConvergence",OpIds::PressureConvergenceEvaluate);
    static const ConservativePressureMethod conservativeMomentum("ConservativePressureMomentum","momentum","rhoU",TargetKind::Physical,
        {OpIds::MomentumAssemble,OpIds::MomentumSolve},{},{"conservativePredictor"});
    static const ConservativePressureMethod conservativeCorrection("ConservativePressureCorrection","pSimple","p",TargetKind::Correction,
        {OpIds::PressureBoundaryPrepare,OpIds::PressureAssemble,OpIds::PressureSolve},{"conservativePredictor"},{"pressureCorrection"});
    static const ConservativePressureMethod conservativeVelocity("ConservativeVelocityCorrection","correctU","rhoU",TargetKind::Physical,
        {OpIds::VelocityCorrect},{"pressureCorrection"});
    static const ConservativePressureMethod conservativeUpdate("ConservativePressureUpdate","correctP","p",TargetKind::Working,
        {OpIds::PressureUpdatePrepare},{"pressureCorrection"},{"preparedPressure"});
    static const ConservativePressureMethod conservativeFlux("ConservativeFluxCorrection","correctFluxp","fluxValidity",TargetKind::Workspace,
        {OpIds::FluxCorrect},{"pressureCorrection"},{"derivedFluxValidity"});
    static const ConservativePressureMethod conservativePublication("ConservativePressurePublication","publishPressure","rhoE",TargetKind::Physical,
        {OpIds::PressureCorrectionCommit},{"preparedPressure","derivedFluxValidity"});
    static const ConservativePressureMethod conservativeRelax("ConservativeFixedTimeRelaxation","relaxIterate","iterate",TargetKind::Workspace,
        {OpIds::PressureRelaxationApply});
    static const ConservativePressureMethod conservativeRestore("ConservativeFluxConsistency","restoreFlux","fluxValidity",TargetKind::Workspace,
        {OpIds::PressureFluxConsistencyRestore},{},{"derivedFluxValidity"});
    static const ConservativePressureMethod conservativeConvergence("ConservativeResidualConvergence","checkConvergence","converged",TargetKind::Workspace,
        {OpIds::PressureConvergenceEvaluate});
    static const TurbulenceTransportMethod turbulence;
    ProviderRegistry result;
    addEulerianMethods(result);
    addImmersedMethods(result);
    result.add(turbulence);
    static const TurbulenceClosureMethod turbulenceClosure;
    result.add(turbulenceClosure);
    result.add(relax);result.add(restore);result.add(convergence);
    for (const auto* provider:{&conservativeMomentum,&conservativeCorrection,&conservativeVelocity,
            &conservativeUpdate,&conservativeFlux,&conservativePublication,&conservativeRelax,
            &conservativeRestore,&conservativeConvergence}) result.add(*provider);
    result.add(residual);
    result.add(evaluation);
    result.add(linear);
    result.add(momentum);
    result.add(correction);
    result.add(pressureUpdate);
    result.add(velocityCorrection);
    result.add(fluxCorrection);
    return result;
}


} // namespace SF::System
