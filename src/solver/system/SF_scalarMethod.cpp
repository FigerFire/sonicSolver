#include "SF_scalarMethod.h"
#include "SF_formulaCompiler.h"
#include "solver/discretization/diffusion/SF_formulaCentral2.h"
#include "core/system/SF_operationIds.h"
#include <cmath>
#include <stdexcept>
namespace SF::System {
namespace {
bool equivalentMath(const FormulaExpr& a,const FormulaExpr& b) {
    return canonicalFormula(Equation{"",a,FormulaExpr::constantValue(0)})==canonicalFormula(Equation{"",b,FormulaExpr::constantValue(0)});
}
class ScalarMethod final : public IProvider {
    bool transport_;
public:
    explicit ScalarMethod(bool transport=false):transport_(transport) {}
    std::string_view id() const override {return transport_?"ScalarTransportUpwindCentral2":"ScalarDiffusionCentral2";}
    std::string_view runtimeProvider() const override {return transport_?ScalarOps::TransportProvider:ScalarOps::Provider;}
    void lowerLifecycle(const ExecutionScope& scope,bool root,
            const std::vector<CompiledEquationCall>&,SolvePlanNode& node) const override {
        if (root && (scope.kind!=ExecutionKind::Sequence || scope.children.empty()
            || scope.children.back().kind!=ExecutionKind::Commit
            || std::count_if(scope.children.begin(),scope.children.end(),[](const auto& item) {
                return item.kind==ExecutionKind::Commit;
            })!=1))
            throw std::runtime_error("Scalar HOW requires one terminal physical Commit.");
        if (scope.kind!=ExecutionKind::Commit) return;
        for (const auto* op:{ScalarOps::Commit,OpIds::TimeCommit}) {
            SolvePlanNode leaf;leaf.kind=PlanNodeKind::Commit;
            leaf.id=leaf.name=leaf.operation=op;leaf.provider=runtimeProvider();
            node.children.push_back(std::move(leaf));
        }
    }
    CompiledEquationCall compile(const ExecutableEquationSystem& system,const EquationCall& call,
            const NumericalBinding& binding) const override {
        if (!binding.inputs.empty() || !binding.parameters.empty())
            throw std::runtime_error("ScalarDiffusionCentral2 has no extra inputs/parameters; D belongs to WHAT.");
        const auto& state=system.state.at(call.target.symbol);
        if (call.target.kind!=TargetKind::Physical || state.shape!=ValueShape::Scalar
            || state.components!=1 || state.componentOffset!=0
            || state.location!=VariableLocation::EulerianCell
            || state.storageBinding!=StorageBinding::NamedDistributed
            || state.constantValue || state.derivation!=StateDerivation::None)
            throw std::runtime_error("Unsupported scalar target: requires writable scalar NamedDistributed cell STATE.");
        const auto& equation=system.registry.at(call.equation);
        using K=FormulaExpr::Kind;
        const auto target=FormulaExpr::symbol(call.target.symbol);
        const auto density=FormulaExpr::symbol("rho");
        const auto concentration=FormulaExpr::divide(target,density);
        auto lhs=FormulaExpr::op("ddt",{target});
        if (transport_) lhs=FormulaExpr::add(lhs,FormulaExpr::op("div",{FormulaExpr::symbol("rhoU"),concentration}));
        if (!equivalentMath(equation.lhs,lhs))
            throw std::runtime_error("Scalar provider LHS does not match its declared transport operator.");
        CompiledScalarEquation contract;contract.target=call.target.symbol;contract.conservativeTransport=transport_;
        contract.stageReads={contract.target};
        if (transport_) {
            const auto& rho=system.state.at("rho");const auto& momentum=system.state.at("rhoU");
            if (rho.components!=1 || momentum.components!=3 || rho.storageBinding!=StorageBinding::PackedDistributed
                || momentum.storageBinding!=StorageBinding::PackedDistributed)
                throw std::runtime_error("Scalar transport requires packed Stage density/momentum.");
            contract.stageReads.push_back("rho");contract.stageReads.push_back("rhoU");
        }
        const auto stageSymbol=[&](const std::string& name) {
            const auto& partner=system.state.at(name);
            if (partner.shape!=ValueShape::Scalar || partner.components!=1 || partner.componentOffset!=0
                || partner.location!=VariableLocation::EulerianCell || partner.constantValue
                || partner.derivation!=StateDerivation::None || partner.storageBinding!=StorageBinding::NamedDistributed)
                throw std::runtime_error("Unsupported scalar Stage partner shape/storage: "+name);
            if (std::find(contract.stageReads.begin(),contract.stageReads.end(),name)==contract.stageReads.end())
                contract.stageReads.push_back(name);
        };
        int diffusion=0;
        const auto validate=[&](const auto& self,const FormulaExpr& node)->void {
            if (node.kind==K::Add && node.arguments.size()==2) {
                self(self,node.arguments[0]);self(self,node.arguments[1]);return;
            }
            if (node.kind==K::Constant && std::isfinite(node.constant)) return;
            if (node.kind==K::Symbol) {stageSymbol(node.name);return;}
            if (node.kind==K::Multiply && node.arguments.size()==2) {
                const auto& a=node.arguments[0];const auto& b=node.arguments[1];
                if (a.kind==K::Constant && std::isfinite(a.constant) && b.kind==K::Symbol) {stageSymbol(b.name);return;}
                if (b.kind==K::Constant && std::isfinite(b.constant) && a.kind==K::Symbol) {stageSymbol(a.name);return;}
                throw std::runtime_error("Unsupported scalar coupling: finite constant times Stage scalar required.");
            }
            if (node.kind!=K::Operator) throw std::runtime_error("Unsupported scalar RHS: only diffusion plus additive known source/constants.");
            if (node.name=="diffusion" && node.arguments.size()==2
                && (transport_?equivalentMath(node.arguments[1],concentration)
                    :node.arguments[1].kind==K::Symbol && node.arguments[1].name==contract.target)) {
                if (++diffusion!=1) throw std::runtime_error("Scalar provider requires one diffusion term.");
                const FormulaExpr* coefficientPtr=&node.arguments[0];
                if (transport_) {
                    const auto& coefficient=*coefficientPtr;
                    if (coefficient.kind!=K::Multiply || coefficient.arguments.size()!=2
                        || !equivalentMath(coefficient.arguments[0],density))
                        throw std::runtime_error("Scalar transport diffusion coefficient requires rho*constantD.");
                    coefficientPtr=&coefficient.arguments[1];
                }
                const auto& coefficient=*coefficientPtr;
                if (coefficient.kind==K::Constant) contract.diffusivity=coefficient.constant;
                else if (coefficient.kind==K::Symbol && system.state.at(coefficient.name).constantValue)
                    contract.diffusivity=*system.state.at(coefficient.name).constantValue;
                else throw std::runtime_error("Unsupported scalar diffusivity: requires constant D.");
                if (!std::isfinite(contract.diffusivity) || contract.diffusivity<=0)
                    throw std::runtime_error("Scalar diffusivity must be finite and positive.");
                return;
            }
            if (node.name=="source" && node.arguments.size()==1 && node.arguments[0].kind==K::Symbol
                && node.arguments[0].name!=contract.target) {
                if (system.state.contains(node.arguments[0].name) && !system.state.at(node.arguments[0].name).constantValue)
                    throw std::runtime_error("Physical partner STATE cannot masquerade as known source callback.");
                contract.sources.push_back(node.arguments[0].name);return;
            }
            throw std::runtime_error("Unsupported scalar operator/operand: "+node.name);
        };
        validate(validate,equation.rhs);
        if (diffusion!=1) throw std::runtime_error("Scalar provider requires one diffusion term.");
        FormulaOperatorCatalog operators;
        operators.add(transport_?Discretization::conservativeScalarCentral2(contract.diffusivity)
            :Discretization::central2ScalarDiffusion(contract.diffusivity));
        FormulaOperatorProvider source;source.id="scalar.known-source";source.mathematicalOperator="source";
        source.matches=[](const auto&) {return true;};
        source.compileValue=[offset=contract.stageReads.size(),names=contract.sources](const FormulaExpr& expression) -> FormulaValueKernel {
            const auto slot=offset+static_cast<std::size_t>(std::find(names.begin(),names.end(),expression.arguments[0].name)-names.begin());
            return [slot](const FormulaValues& values,int cell,int) {
                return values.boundReads.at(slot)(cell,0);
            };
        };
        operators.add(std::move(source));
        CompiledEquationCall result;
        contract.rhs=compileFormulaValue(equation,equation.rhs,operators,result.operatorBindings,contract.stageReads);
        if(transport_)result.operatorBindings.push_back("scalar.conservative.upwind1");
        const auto stageReads=contract.stageReads;
        result.providerContract=std::move(contract);
        result.source=call;result.target.symbol=call.target.symbol;result.target.kind=call.target.kind;
        result.equationMethod=id();result.temporalResidual=true;result.synchronousStages=true;
        result.temporalStepSize=ScalarOps::Dt;
        result.temporalSnapshot=instanceOperation(ScalarOps::Begin,call.occurrence);
        result.temporalStagePrepare=instanceOperation("scalar.stage.prepare",call.occurrence);
        result.temporalRhs=instanceOperation("scalar.stage.rhs",call.occurrence);
        result.temporalAdvance=instanceOperation("scalar.stage.advance",call.occurrence);
        result.temporalPublish=instanceOperation(ScalarOps::Commit,call.occurrence);
        result.backendOperation=ScalarOps::Stage;
        result.oldTimeWorkspace="scalar."+call.target.symbol+".old";
        result.stageWorkspace="scalar."+call.target.symbol+".stage";
        result.residualWorkspace="scalar."+call.target.symbol+".rhs";
        result.requiredHaloWidth=1;
        result.numericalWorkspaces={result.oldTimeWorkspace,result.stageWorkspace,result.residualWorkspace};
        result.temporalCapabilities={FDM::TimeRecipeId::ForwardEuler,FDM::TimeRecipeId::ClassicalRK4};
        result.calls={{call.equation,call.target.symbol}};
        result.reads=stageReads;result.writes={call.target.symbol};
        for (const auto& name:stageReads) result.stateUses.push_back({name,StateVersion::Stage});
        result.stateEffects={{result.residualWorkspace,true,false,false}};
        result.requirements={"serial single patch orthogonal grid","common group dt","scalar fixedValue/zeroGradient boundary","stage-time known source"};
        if(transport_)result.requirements.push_back("conservative Upwind1 advection and rho-weighted Central2 concentration diffusion; boundary values refer to rhoC");
        for (const auto& op:{result.temporalSnapshot,result.temporalStagePrepare,result.temporalRhs,
                result.temporalAdvance,result.temporalPublish})
            result.operations.push_back({op,op,op==result.temporalPublish?OperationStage::StepCommit:OperationStage::Prepare,
                {OperationCapability::ScalarExplicit},{OriginKind::Generated,std::string(id())}});
        return result;
    }
};
}
void addScalarMethod(ProviderRegistry& registry) {static const ScalarMethod method,transport(true);registry.add(method);registry.add(transport);}
}
