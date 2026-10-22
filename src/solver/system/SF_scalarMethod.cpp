#include "SF_builtinState.h"
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
    bool shared_;
    bool model_;
public:
    explicit ScalarMethod(bool transport=false,bool shared=false,bool model=false):transport_(transport),shared_(shared),model_(model) {}
    std::string_view id() const override {return model_?"SSTStageTransport":shared_?"ScalarTransportSharedFluxCentral2":transport_?"ScalarTransportUpwindCentral2":"ScalarDiffusionCentral2";}
    std::string_view runtimeProvider() const override {return model_?ScalarOps::StageTransportProvider:shared_?ScalarOps::SharedFluxProvider:transport_?ScalarOps::TransportProvider:ScalarOps::Provider;}
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
        if ((!model_ && !binding.inputs.empty()) || !binding.parameters.empty())
            throw std::runtime_error("ScalarDiffusionCentral2 has no extra inputs/parameters; D belongs to WHAT.");
        const auto& state=system.state.at(call.target.symbol);
        if (call.target.kind!=TargetKind::Physical || state.shape!=ValueShape::Scalar
            || state.components!=1 || state.componentOffset!=0
            || state.location!=VariableLocation::EulerianCell
            || (model_?state.storageBinding!=StorageBinding::ProviderDistributed:state.storageBinding!=StorageBinding::NamedDistributed)
            || state.constantValue || state.readOnly())
            throw std::runtime_error("Unsupported scalar target: requires writable scalar NamedDistributed cell STATE.");
        const auto& equation=system.registry.at(call.equation);
        using K=FormulaExpr::Kind;
        const auto target=state.expression();
        const auto density=BuiltinStateCatalog().expression("rho");
        const auto concentration=model_?target:FormulaExpr::divide(target,density);
        auto lhs=FormulaExpr::op(Math::Op::Ddt,{model_?FormulaExpr::multiply(density,target):target});
        if (transport_) lhs=FormulaExpr::add(lhs,FormulaExpr::op(Math::Op::Div,{BuiltinStateCatalog().expression("rhoU"),concentration}));
        if (transport_ && !equivalentMath(equation.lhs,lhs))
            throw std::runtime_error("Scalar provider LHS does not match its declared transport operator.");
        CompiledScalarEquation contract;contract.target=call.target.symbol;contract.conservativeTransport=transport_;contract.sharedMassFlux=shared_;contract.primitiveStorage=contract.modelCoefficients=model_;
        if(model_) {
            if(binding.inputs!=std::vector<std::string>{"k","omega"} || (call.target.symbol!="k" && call.target.symbol!="omega"))throw std::runtime_error("Stage SST requires explicit k/omega coefficient port.");
            contract.modelTargets=binding.inputs;
        }
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
                || partner.readOnly() || (model_?partner.storageBinding!=StorageBinding::ProviderDistributed && partner.storageBinding!=StorageBinding::TransientWorkspace:partner.storageBinding!=StorageBinding::NamedDistributed))
                throw std::runtime_error("Unsupported scalar Stage partner shape/storage: "+name);
            if (std::find(contract.stageReads.begin(),contract.stageReads.end(),name)==contract.stageReads.end())
                contract.stageReads.push_back(name);
        };
        const auto rhs=transport_?equation.rhs:explicitFormulaRhs(equation,contract.target);
        const auto structure=inspectFormula(rhs,contract.target);
        int diffusion=0;
        for(const auto& item:structure.nodes) {
            const auto& node=*item.expression;
            if (node.kind()==K::Operator && node.operation()==Math::Op::Diffusion) {
                if(node.arguments().size()!=2 || !(transport_?equivalentMath(node.arguments()[1],concentration)
                    :node.arguments()[1].kind()==K::Symbol && node.arguments()[1].name()==contract.target))
                    throw std::runtime_error("Unsupported: selected Central2 kernel only diffuses its declared target/concentration");
                if (++diffusion!=1) throw std::runtime_error("Scalar provider requires one diffusion term.");
                const FormulaExpr* coefficientPtr=&node.arguments()[0];
                if (transport_) {
                    const auto& coefficient=*coefficientPtr;
                    if (coefficient.kind()!=K::Multiply || coefficient.arguments().size()!=2
                        || !equivalentMath(coefficient.arguments()[0],density))
                        throw std::runtime_error("Scalar transport diffusion coefficient requires rho*constantD.");
                    coefficientPtr=&coefficient.arguments()[1];
                }
                const auto& coefficient=*coefficientPtr;
                if (coefficient.kind()==K::Constant) contract.diffusivity=coefficient.constant();
                else if (coefficient.kind()==K::Symbol && system.state.at(coefficient.name()).constantValue)
                    contract.diffusivity=*system.state.at(coefficient.name()).constantValue;
                else if(shared_ && coefficient.kind()==K::Symbol) {
                    stageSymbol(coefficient.name());contract.diffusivityState=coefficient.name();
                } else throw std::runtime_error("Unsupported scalar diffusivity: requires constant D or explicit shared-flux Stage coefficient.");
                if (contract.diffusivityState.empty() && (!std::isfinite(contract.diffusivity) || (shared_?contract.diffusivity<0:contract.diffusivity<=0)))
                    throw std::runtime_error("Scalar diffusivity must be finite and positive.");
                continue;
            }
            if (node.kind()==K::Operator && node.operation()==Math::Op::Source) {
                if(node.arguments().size()!=1 || node.arguments()[0].kind()!=K::Symbol || node.arguments()[0].name()==contract.target)
                    throw std::runtime_error("Unsupported: known-source provider requires one declared callback symbol");
                if (system.state.contains(node.arguments()[0].name()) && !system.state.at(node.arguments()[0].name()).constantValue)
                    throw std::runtime_error("Physical partner STATE cannot masquerade as known source callback.");
                if(std::find(contract.sources.begin(),contract.sources.end(),node.arguments()[0].name())==contract.sources.end())contract.sources.push_back(node.arguments()[0].name());
            }
        }
        for(const auto& name:structure.symbols) {
            if(std::find(contract.stageReads.begin(),contract.stageReads.end(),name)!=contract.stageReads.end())continue;
            if(std::find(contract.sources.begin(),contract.sources.end(),name)!=contract.sources.end())continue;
            if(system.state.contains(name) && system.state.at(name).constantValue)continue;
            stageSymbol(name);
        }
        if (diffusion!=1) throw std::runtime_error("Scalar provider requires one diffusion term.");
        FormulaOperatorCatalog operators;
        if(shared_) {
            auto slot=std::find(contract.stageReads.begin(),contract.stageReads.end(),contract.diffusivityState)-contract.stageReads.begin();
            operators.add(Discretization::conservativeScalarCentral2(contract.diffusivity,contract.diffusivityState.empty()?-1:static_cast<int>(slot),model_));
        } else operators.add(transport_?Discretization::conservativeScalarCentral2(contract.diffusivity)
            :Discretization::central2ScalarDiffusion(contract.diffusivity));
        FormulaOperatorProvider source;source.id="scalar.known-source";source.mathematicalOperator=Math::Op::Source;
        source.matches=[](const auto&) {return true;};
        source.compileValue=[offset=contract.stageReads.size(),names=contract.sources](const FormulaExpr& expression) -> FormulaValueKernel {
            const auto slot=offset+static_cast<std::size_t>(std::find(names.begin(),names.end(),expression.arguments()[0].name())-names.begin());
            return [slot](const FormulaValues& values,int cell,int) {
                return values.boundReads.at(slot)(cell,0);
            };
        };
        operators.add(std::move(source));
        CompiledEquationCall result;
        contract.nonlinear=structure.targetDegree>1;
        contract.normalized=formulaText(transport_?equation:Equation{equation.id,FormulaExpr::op(Math::Op::Ddt,{target}),rhs,equation.origin});
        contract.rhs=compileFormulaValue(equation,rhs,operators,result.operatorBindings,contract.stageReads,&system.state,&contract.lowering);
        if(transport_)result.operatorBindings.push_back(shared_?"scalar.canonical-mass-flux.upwind1":"scalar.independent-center-rhoU.upwind1");
        if(shared_)result.stageRequires={"canonicalMassFlux(Stage)"};
        if(model_) {
            result.stageRequires.push_back("transportStageCoefficients");
            if(call.target.symbol=="k") {
                result.coefficientPrepare=instanceOperation("scalar.stage.coefficients",call.occurrence);
                result.coefficientPublish=instanceOperation("scalar.coefficients.publish",call.occurrence);
                result.coefficientProvides={"transportStageCoefficients"};
                result.coefficientStates={"mu_t","SST.nu.k","SST.nu.omega","SST.production.k","SST.production.omega","SST.sink.k","SST.sink.omega","SST.cross"};
                for(const auto& name:result.coefficientStates)system.state.at(name);
            }
        }
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
        if(model_)for(const auto& name:stageReads)if(name.rfind("SST.",0)==0)
            result.coefficientUses.push_back({name,{StateVersion::Stage}});
        result.stateEffects={{result.residualWorkspace,true,false,false}};
        result.requirements={"serial single patch orthogonal grid","common group dt","scalar fixedValue/zeroGradient boundary","stage-time known source"};
        if(model_)result.requirements.push_back("primitive physical alias; conservative stage workspace; explicit model coefficient publication; no correct() advance");
        if(transport_)result.requirements.push_back("conservative Upwind1 advection and rho-weighted Central2 concentration diffusion; boundary values refer to rhoC");
        for (const auto& op:{result.temporalSnapshot,result.temporalStagePrepare,result.temporalRhs,
                result.temporalAdvance,result.temporalPublish})
            result.operations.push_back({op,op,op==result.temporalPublish?OperationStage::StepCommit:OperationStage::Prepare,
                {OperationCapability::ScalarExplicit},{OriginKind::Generated,std::string(id())}});
        for(const auto& op:{result.coefficientPrepare,result.coefficientPublish})if(!op.empty())
            result.operations.push_back({op,op,op==result.coefficientPublish?OperationStage::StepCommit:OperationStage::Prepare,
                {OperationCapability::ScalarExplicit},{OriginKind::Generated,std::string(id())}});
        return result;
    }
};
}
void addScalarMethod(ProviderRegistry& registry) {static const ScalarMethod method,transport(true),shared(true,true);registry.add(method);registry.add(transport);registry.add(shared);static const ScalarMethod stage(true,true,true);registry.add(stage);}
}
