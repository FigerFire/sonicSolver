#include "core/interfaces/SF_scalarEquationData.h"
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
ResourceViewRequirement fieldResource(std::string id,const std::string& target,
        const std::vector<StateUse>& reads,std::vector<std::string> boundaries={},
        std::vector<std::string> sources={},bool local=false) {
    FieldDataRequirement fields;fields.target=target;fields.states={target};fields.localOnly=local;
    for(const auto& read:reads)fields.states.push_back(read.symbol);
    fields.boundaries=std::move(boundaries);fields.sources=std::move(sources);
    return {std::move(id),"native.scalar",typeid(FDM::ScalarEquationData),false,std::move(fields)};
}
enum class ScalarTransportForm { None, Conservative, SharedConservative, SharedPrimitive };
class ScalarMethod final : public IProvider {
    ScalarTransportForm form_;
    bool arithmetic_;
public:
    explicit ScalarMethod(ScalarTransportForm form=ScalarTransportForm::None,bool arithmetic=false):form_(form),arithmetic_(arithmetic) {}
    std::string_view id() const override {
        return form_==ScalarTransportForm::SharedPrimitive?"CapacityTransportSharedFluxCentral2"
            :form_==ScalarTransportForm::SharedConservative?"ScalarTransportSharedFluxCentral2"
            :form_==ScalarTransportForm::Conservative?"ScalarTransportUpwindCentral2"
            :arithmetic_?"ExplicitCentral2Arithmetic":"ScalarDiffusionCentral2";
    }
    std::string_view runtimeProvider() const override {
        return form_==ScalarTransportForm::SharedPrimitive?ScalarOps::CapacityProvider
            :form_==ScalarTransportForm::SharedConservative?ScalarOps::SharedFluxProvider
            :form_==ScalarTransportForm::Conservative?ScalarOps::TransportProvider:ScalarOps::Provider;
    }
    CompiledEquationCall compile(const ExecutableEquationSystem& system,const EquationCall& call,
            const NumericalBinding& binding) const override {
        const bool transport=form_!=ScalarTransportForm::None;
        const bool shared=form_==ScalarTransportForm::SharedConservative || form_==ScalarTransportForm::SharedPrimitive;
        const bool primitive=form_==ScalarTransportForm::SharedPrimitive;
        if (!binding.inputs.empty() || !binding.parameters.empty())
            throw std::runtime_error("ScalarDiffusionCentral2 has no extra inputs/parameters; D belongs to WHAT.");
        const auto& state=system.state.at(call.target.symbol);
        if (call.target.kind!=TargetKind::Physical || (state.shape!=ValueShape::Scalar && (transport || state.shape!=ValueShape::Vector || state.components!=3))
            || (state.shape==ValueShape::Scalar && state.components!=1) || state.componentOffset!=0
            || state.location!=VariableLocation::EulerianCell
            || (state.storageBinding!=StorageBinding::NamedDistributed && (!primitive || state.storageBinding!=StorageBinding::ProviderDistributed))
            || state.constantValue || state.readOnly())
            throw std::runtime_error("Unsupported scalar target: requires writable scalar or three-component vector NamedDistributed cell STATE.");
        const auto& equation=system.registry.at(call.equation);
        using K=FormulaExpr::Kind;
        const auto target=state.expression();
        auto concentration=target;
        if(transport) {
            const auto density=BuiltinStateCatalog().expression("rho");
            concentration=primitive?target:FormulaExpr::divide(target,density);
            const auto transient=FormulaExpr::op(Math::Op::Ddt,{primitive?FormulaExpr::multiply(density,target):target});
            const auto lhs=FormulaExpr::add(transient,FormulaExpr::op(Math::Op::Div,{BuiltinStateCatalog().expression("rhoU"),concentration}));
            if(!equivalentMath(equation.lhs,lhs))
                throw std::runtime_error("Selected conservative transport implementation requires its declared capacity and mass-flux LHS; generic explicit mass inversion is not implemented.");
        }
        double transportD=0;
        CompiledScalarEquation contract;contract.components=state.components;contract.arithmeticFaces=arithmetic_;contract.target=call.target.symbol;contract.conservativeTransport=transport;contract.sharedMassFlux=shared;contract.primitiveStorage=primitive;
        contract.valueReads={{contract.target,StateVersion::Stage}};
        if (transport) {
            const auto& rho=system.state.at("rho");const auto& momentum=system.state.at("rhoU");
            if (rho.components!=1 || momentum.components!=3 || rho.storageBinding!=StorageBinding::PackedDistributed
                || momentum.storageBinding!=StorageBinding::PackedDistributed)
                throw std::runtime_error("Scalar transport requires packed Stage density/momentum.");
            contract.valueReads.push_back({"rho",StateVersion::Stage});contract.valueReads.push_back({"rhoU",StateVersion::Stage});
        }
        const auto appendRead=[&](const std::string& name,StateVersion version) {
            const auto found=std::find_if(contract.valueReads.begin(),contract.valueReads.end(),[&](const auto& r){return r.symbol==name;});
            if(found!=contract.valueReads.end()) {
                if(name==contract.target && version==StateVersion::OldTime && !transport) {found->version=version;return;}
                if(found->version!=version)throw std::runtime_error("Conflicting value read versions: "+name);
            } else contract.valueReads.push_back({name,version});
        };
        const auto readSymbol=[&](const auto& self,const std::string& name,std::optional<StateVersion> inherited=std::nullopt)->void {
            const auto& partner=system.state.at(name);
            const auto explicitRead=std::find_if(call.inputs.begin(),call.inputs.end(),[&](const auto& r){return r.symbol==name;});
            const bool evolved=std::find(system.state.solutionVariables().begin(),system.state.solutionVariables().end(),name)!=system.state.solutionVariables().end();
            if(!transport && !inherited && !evolved && explicitRead==call.inputs.end())throw std::runtime_error("Missing explicit input version: "+name);
            const auto version=inherited.value_or(explicitRead==call.inputs.end()?StateVersion::Stage:explicitRead->version);
            if(inherited && explicitRead!=call.inputs.end() && explicitRead->version!=version)throw std::runtime_error("Producer dependency version conflict: "+name);
            if(version!=StateVersion::Current && version!=StateVersion::Stage && version!=StateVersion::OldTime)throw std::runtime_error("Unsupported input version: "+name+" / "+toString(version)+" requires a real snapshot binding");
            if(version==StateVersion::OldTime && (transport || primitive || !evolved || partner.producer))throw std::runtime_error("Unsupported OldTime input: requires an existing ordinary temporal snapshot: "+name);
            // The complete program validates actual Stage publishers; membership in
            // the evolved unknown set does not decide whether a coefficient is readable.
            if((partner.shape!=ValueShape::Scalar && (transport || partner.shape!=ValueShape::Vector || partner.components!=3)) || partner.componentOffset!=0
                || partner.location!=VariableLocation::EulerianCell || partner.constantValue)
                throw std::runtime_error("Unsupported scalar Stage partner shape/storage: "+name);
            if(partner.producer) {
                if(!arithmetic_ || partner.evaluation!=StateEvaluation::Lazy)throw std::runtime_error("Unsupported coefficient producer: explicit materialized HOW equation/publication required: "+name);
                const bool visited=std::any_of(contract.valueReads.begin(),contract.valueReads.end(),[&](const auto& r){return r.symbol==name;});
                appendRead(name,version);if(visited)return;
                for(const auto& dep:partner.dependencies)self(self,dep,version);
                return;
            }
            if(partner.storageBinding!=StorageBinding::NamedDistributed && (!(primitive || arithmetic_) || (partner.storageBinding!=StorageBinding::ProviderDistributed && partner.storageBinding!=StorageBinding::TransientWorkspace)))
                throw std::runtime_error("Unsupported scalar Stage partner shape/storage: "+name);
            appendRead(name,version);
        };
        const auto stageSymbol=[&](const std::string& name){readSymbol(readSymbol,name);};
        const auto rhs=transport?equation.rhs:explicitFormulaRhs(equation,contract.target);
        const auto structure=inspectFormula(rhs,contract.target);
        int diffusion=0;
        for(const auto& item:structure.nodes) {
            const auto& node=*item.expression;
            if (node.kind()==K::Operator && node.operation()==Math::Op::Diffusion) {
                if(node.arguments().size()!=2 || !(transport?equivalentMath(node.arguments()[1],concentration)
                    :node.arguments()[1].kind()==K::Symbol))
                    throw std::runtime_error("Unsupported: Central2 diffusion requires a declared Stage symbol/concentration operand");
                if (++diffusion!=1 && transport) throw std::runtime_error("Scalar transport provider requires one diffusion term.");
                if(!transport)stageSymbol(node.arguments()[1].name());
                std::string coefficientState;
                const FormulaExpr* coefficientPtr=&node.arguments()[0];
                if (transport) {
                    const auto& coefficient=*coefficientPtr;
                    if (coefficient.kind()!=K::Multiply || coefficient.arguments().size()!=2
                        || !equivalentMath(coefficient.arguments()[0],BuiltinStateCatalog().expression("rho")))
                        throw std::runtime_error("Scalar transport diffusion coefficient requires rho*constantD.");
                    coefficientPtr=&coefficient.arguments()[1];
                }
                const auto& coefficient=*coefficientPtr;
                if (coefficient.kind()==K::Constant) transportD=coefficient.constant();
                else if (coefficient.kind()==K::Symbol && system.state.at(coefficient.name()).constantValue)
                    transportD=*system.state.at(coefficient.name()).constantValue;
                else if(arithmetic_ && coefficient.kind()==K::Symbol) {
                    stageSymbol(coefficient.name());coefficientState=coefficient.name();
                }
                else if(shared && coefficient.kind()==K::Symbol) {
                    stageSymbol(coefficient.name());contract.diffusivityState=coefficient.name();
                    const auto& field=system.state.at(coefficient.name());
                    if(field.storageBinding==StorageBinding::TransientWorkspace && field.workspaceBoundary!=WorkspaceBoundary::NearestInterior)
                        throw std::runtime_error("Missing coefficient read-region preparation: step="+call.occurrence+" equation="+call.equation+" target="+call.target.symbol+" input="+field.id);
                } else throw std::runtime_error("Unsupported scalar diffusivity: requires constant D or explicit shared-flux Stage coefficient.");
                if (contract.diffusivityState.empty() && (!std::isfinite(transportD) || (transport && (shared?transportD<0:transportD<=0))))
                    throw std::runtime_error("Diffusion coefficient must be finite; transport requires positive/nonnegative diffusivity.");
                contract.diffusions.push_back({node.occurrence().empty()?item.path:node.occurrence(),
                    transport?contract.target:node.arguments()[1].name(),transportD,node.arguments()[1].signature().rank==0?1:3,coefficientState});
                continue;
            }
            if (node.kind()==K::Operator && node.operation()==Math::Op::Source) {
                if(contract.components!=1)throw std::runtime_error("Unsupported: vector known-source callbacks need a component-aware source port.");
                if(node.arguments().size()!=1 || node.arguments()[0].kind()!=K::Symbol || node.arguments()[0].name()==contract.target)
                    throw std::runtime_error("Unsupported: known-source provider requires one declared callback symbol");
                if (system.state.contains(node.arguments()[0].name()) && !system.state.at(node.arguments()[0].name()).constantValue)
                    throw std::runtime_error("Physical partner STATE cannot masquerade as known source callback.");
                if(std::find(contract.sources.begin(),contract.sources.end(),node.arguments()[0].name())==contract.sources.end())contract.sources.push_back(node.arguments()[0].name());
            }
        }
        for(const auto& name:structure.symbols) {
            if(std::any_of(contract.valueReads.begin(),contract.valueReads.end(),[&](const auto& r){return r.symbol==name;}))continue;
            if(std::find(contract.sources.begin(),contract.sources.end(),name)!=contract.sources.end())continue;
            if(system.state.contains(name) && system.state.at(name).constantValue)continue;
            stageSymbol(name);
        }
        if (transport && diffusion!=1) throw std::runtime_error("Scalar transport provider requires one diffusion term.");
        if(arithmetic_)for(const auto& d:contract.diffusions) {
            const auto input=std::find_if(contract.valueReads.begin(),contract.valueReads.end(),[&](const auto& r){return r.symbol==d.input;});
            if(input==contract.valueReads.end() || input->version!=StateVersion::Stage)
                throw std::runtime_error("Unsupported diffusion operand version: published Stage unknown required: "+d.input);
        }
        for(const auto& read:call.inputs) {
            if(read.kind!=TargetKind::Physical || read.halo || !read.snapshot.empty() || read.perIteration)throw std::runtime_error("Unsupported authored scalar read qualifier: "+read.symbol);
            if(std::count_if(call.inputs.begin(),call.inputs.end(),[&](const auto& r){return r.symbol==read.symbol;})!=1)throw std::runtime_error("Duplicate input read declaration: "+read.symbol);
            if(std::none_of(contract.valueReads.begin(),contract.valueReads.end(),[&](const auto& r){return r.symbol==read.symbol && r.version==read.version;}))throw std::runtime_error("Unconsumed input version declaration: "+read.symbol);
        }
        std::vector<std::string> slots;for(const auto& read:contract.valueReads)slots.push_back(read.symbol);
        FormulaOperatorCatalog operators;
        if(shared) {
            auto slot=std::find(slots.begin(),slots.end(),contract.diffusivityState)-slots.begin();
            operators.add(Discretization::conservativeScalarCentral2(transportD,contract.diffusivityState.empty()?-1:static_cast<int>(slot),primitive));
        } else if(transport)operators.add(Discretization::conservativeScalarCentral2(transportD));
        else operators.add(Discretization::central2ScalarDiffusion(system.state,slots,arithmetic_));
        FormulaOperatorProvider source;source.id="scalar.known-source";source.mathematicalOperator=Math::Op::Source;
        source.matches=[](const auto&) {return true;};
        source.compileValue=[offset=slots.size(),names=contract.sources](const FormulaExpr& expression) -> FormulaValueKernel {
            const auto slot=offset+static_cast<std::size_t>(std::find(names.begin(),names.end(),expression.arguments()[0].name())-names.begin());
            return [slot](const FormulaValues& values,int cell,int) {
                return values.boundReads.at(slot)(cell,0);
            };
        };
        operators.add(std::move(source));
        CompiledEquationCall result;
        contract.nonlinear=structure.targetDegree>1;
        contract.normalized=formulaText(transport?equation:Equation{equation.id,FormulaExpr::op(Math::Op::Ddt,{target}),rhs,equation.origin});
        contract.rhs=compileFormulaValue(equation,rhs,operators,result.operatorBindings,slots,&system.state,&contract.lowering,binding.operators);
        if(transport)result.operatorBindings.push_back(shared?"scalar.canonical-mass-flux.upwind1":"scalar.independent-center-rhoU.upwind1");
        if(shared)result.stageRequires={"canonicalMassFlux(Stage)"};
        const auto valueReads=contract.valueReads;
        {
            std::vector<std::string> boundaries;
            if(diffusion || transport)boundaries.push_back(call.target.symbol);
            for(const auto& d:contract.diffusions)boundaries.push_back(d.input);
            auto resource=fieldResource(call.occurrence,call.target.symbol,contract.valueReads,std::move(boundaries),contract.sources);
            if(state.storageBinding==StorageBinding::ProviderDistributed)resource.owner="model.scalar";
            result.physicalResources.push_back(std::move(resource));
        }
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
        result.requiredHaloWidth=diffusion?1:0;
        result.numericalWorkspaces={result.oldTimeWorkspace,result.stageWorkspace,result.residualWorkspace};
        result.temporalCapabilities={FDM::TimeRecipeId::ForwardEuler,FDM::TimeRecipeId::ClassicalRK4};
        result.calls={{call.equation,call.target.symbol}};
        for(const auto& read:valueReads)result.reads.push_back(read.symbol);result.writes={call.target.symbol};
        result.stateUses=valueReads;
        result.stateEffects={{result.residualWorkspace,true,false,false}};
        result.requirements={"serial single patch","common group dt","component-major scalar/vector workspace","stage-time known source"};
        if(diffusion)result.requirements.push_back(arithmetic_
            ? "Cartesian Central2 arithmetic faces; existing nodal fixedValue/zeroGradient ghost; coefficient neighborhood reads mirrored physical nodes in the selected version"
            : "Cartesian Central2; each stencil input requires its own fixedValue/zeroGradient Stage boundary");
        if(primitive)result.requirements.push_back("primitive physical alias; rho-weighted conservative stage workspace; declared coefficient publishers");
        if(transport)result.requirements.push_back("conservative Upwind1 advection and rho-weighted Central2 concentration diffusion; boundary values refer to the declared physical scalar target");
        for (const auto& op:{result.temporalSnapshot,result.temporalStagePrepare,result.temporalRhs,
                result.temporalAdvance,result.temporalPublish})
            result.operations.push_back({op,op,op==result.temporalPublish?OperationStage::StepCommit:OperationStage::Prepare,
                {OperationCapability::ScalarExplicit},{OriginKind::Generated,std::string(id())}});
        return result;
    }
};
/// Exact local affine decomposition; never numerical differencing at target=0/1.
 } // namespace
CompiledEquationCall compileLocalAffine(const ExecutableEquationSystem& system,const EquationCall& call,const NumericalBinding& binding) {
        if(!binding.inputs.empty() || !binding.parameters.empty() || !binding.operators.empty())throw std::runtime_error("LocalAffine accepts STATE read versions, not hidden method/operator parameters");
        const auto& target=system.state.at(call.target.symbol);
        if(call.target.kind==TargetKind::Workspace || (call.stageScoped && call.target.kind!=TargetKind::Physical) || target.readOnly() || target.shape!=ValueShape::Scalar || target.components!=1 || target.componentOffset!=0 || target.location!=VariableLocation::EulerianCell || (target.storageBinding!=StorageBinding::NamedDistributed && target.storageBinding!=StorageBinding::ProviderDistributed && target.storageBinding!=StorageBinding::TransientWorkspace))throw std::runtime_error("LocalAffine requires a writable scalar cell target with an existing or solver workspace binding");
        const auto& equation=system.registry.at(call.equation);const auto residual=FormulaExpr::subtract(equation.lhs,equation.rhs);const auto shape=inspectFormula(residual,call.target.symbol);
        if(std::find(shape.symbols.begin(),shape.symbols.end(),call.target.symbol)==shape.symbols.end())throw std::runtime_error("HOW target is absent from equation: "+call.target.symbol);
        if(shape.targetDegree!=1)throw std::runtime_error("Unsupported LocalAffine: nonlinear or target-dependent denominator; select an implemented nonlinear method with branch/convergence control");
        for(auto operation:shape.operations)if(operation!=Math::Op::Add && operation!=Math::Op::Subtract && operation!=Math::Op::Multiply && operation!=Math::Op::Divide && operation!=Math::Op::Negate && operation!=Math::Op::Pow && operation!=Math::Op::Component && operation!=Math::Op::Vector)throw std::runtime_error("LocalAffine supports local algebra only; no differential/operator solve is implemented");
        CompiledAffineEquation contract;contract.target=target.id;contract.boundaryPreparation=target.workspaceBoundary;std::vector<std::string> slots;
        for(const auto& name:shape.symbols)if(name!=target.id && !system.state.at(name).constantValue) {
            const auto& input=system.state.at(name);
            if((input.components!=1 && !(input.components==3 && input.shape==ValueShape::Vector)) || input.location!=VariableLocation::EulerianCell || input.componentOffset<0
                || (input.storageBinding!=StorageBinding::NamedDistributed && input.storageBinding!=StorageBinding::ProviderDistributed
                    && input.storageBinding!=StorageBinding::PackedDistributed && input.storageBinding!=StorageBinding::TransientWorkspace) || input.producer)
                throw std::runtime_error("Unsupported LocalAffine input resource/layout: "+name);
            StateUse use{name,call.stageScoped?StateVersion::Stage:StateVersion::Current};
            const auto found=std::find_if(call.inputs.begin(),call.inputs.end(),[&](const auto& r){return r.symbol==name;});
            if(found!=call.inputs.end())use=*found;
            if((use.version!=StateVersion::Current && use.version!=StateVersion::Stage && use.version!=StateVersion::OldTime) || use.halo
                || !use.snapshot.empty() || use.perIteration || (!call.stageScoped && use.version!=StateVersion::Current)
                || use.kind==TargetKind::Workspace || (use.kind!=TargetKind::Physical && (call.stageScoped || use.version!=StateVersion::Current)))
                throw std::runtime_error("Unsupported LocalAffine input version/snapshot: "+name);
            contract.valueReads.push_back(use);slots.push_back(name);
        }
        for(const auto& read:call.inputs)
            if(std::count_if(call.inputs.begin(),call.inputs.end(),[&](const auto& r){return r.symbol==read.symbol;})!=1
                || std::none_of(contract.valueReads.begin(),contract.valueReads.end(),[&](const auto& r){return r.symbol==read.symbol;}))
                throw std::runtime_error("Duplicate/unconsumed LocalAffine input (target remains an unknown): "+read.symbol);
        using Kernel=FormulaValueKernel;struct Pair{Kernel a,r;};const auto constant=[](double x)->Kernel{return [x](const FormulaValues&,int,int){return x;};};
        const auto combine=[](Kernel left,Kernel right,char op)->Kernel{return [left=std::move(left),right=std::move(right),op](const FormulaValues& v,int cell,int c){double a=left(v,cell,c),b=right(v,cell,c);if(op=='+')return a+b;if(op=='-')return a-b;if(op=='*')return a*b;if(b==0)throw std::runtime_error("LocalAffine coefficient expression divides by zero");return a/b;};};
        const auto split=[&](const auto& self,const FormulaExpr& e)->Pair {
            const auto part=inspectFormula(e,target.id);if(part.targetDegree==0){std::vector<std::string> used;return {constant(0),compileFormulaValue(equation,e,{},used,slots,&system.state)};}
            if(e.kind()==FormulaExpr::Kind::Symbol)return {constant(1),constant(0)};
            if(e.kind()==FormulaExpr::Kind::Negate){auto p=self(self,e.arguments()[0]);return {combine(constant(-1),p.a,'*'),combine(constant(-1),p.r,'*')};}
            if(e.operation()==Math::Op::Pow && e.arguments()[1].constant()==1)return self(self,e.arguments()[0]);
            auto left=self(self,e.arguments()[0]),right=self(self,e.arguments()[1]);
            if(e.kind()==FormulaExpr::Kind::Add || e.kind()==FormulaExpr::Kind::Subtract){char op=e.kind()==FormulaExpr::Kind::Add?'+':'-';return {combine(left.a,right.a,op),combine(left.r,right.r,op)};}
            if(e.kind()==FormulaExpr::Kind::Divide)return {combine(left.a,right.r,'/'),combine(left.r,right.r,'/')};
            if(e.kind()==FormulaExpr::Kind::Multiply){bool l=inspectFormula(e.arguments()[0],target.id).targetDegree!=0;return {combine(l?left.a:left.r,l?right.r:right.a,'*'),combine(left.r,right.r,'*')};}
            throw std::runtime_error("Unsupported LocalAffine target structure");
        };
        auto pair=split(split,residual);contract.coefficient=std::move(pair.a);contract.remainder=std::move(pair.r);
        CompiledEquationCall result;result.source=call;result.target.symbol=target.id;result.target.kind=call.target.kind;result.stageLocal=call.stageScoped;result.equationMethod="LocalAffine";result.backendProvider=ScalarOps::Provider;result.calls={{call.equation,target.id}};result.writes={target.id};result.stateUses=contract.valueReads;for(const auto& input:contract.valueReads)result.reads.push_back(input.symbol);result.stateEffects={{target.id,true,false,true}};
        if(call.target.kind!=TargetKind::Physical) {
            result.target.workspace="state."+target.id+(call.target.kind==TargetKind::Working?".Working":".Correction");
            result.stateEffects.clear();result.workspaceProvides={result.target.workspace};
            result.stepWorkspaces=result.workspaceProvides;
        }
        for(const auto& input:contract.valueReads)if(input.kind!=TargetKind::Physical)
            result.workspaceRequires.push_back("state."+input.symbol+(input.kind==TargetKind::Working?".Working":".Correction"));
        result.residualCapabilities={ResidualKind::StateChange,ResidualKind::DiscreteEquation};
        result.backendOperation=instanceOperation("affine.solve",call.occurrence);result.operations.push_back({result.backendOperation,"solve local affine residual for "+target.id,OperationStage::Prepare,{OperationCapability::ScalarExplicit},{OriginKind::Generated,"LocalAffine"}});
        if(result.stageLocal){result.stageWorkspace="affine."+call.occurrence+".stage";result.numericalWorkspaces={result.stageWorkspace};result.temporalPublish=instanceOperation("affine.publish",call.occurrence);result.operations.push_back({result.temporalPublish,"publish algebraic relation from committed physical inputs",OperationStage::StepCommit,{OperationCapability::ScalarExplicit},{OriginKind::Generated,"LocalAffine"}});}
        result.requirements={"serial single patch","local scalar affine residual; no root branch selection",result.stageLocal?"explicit HOW Stage frequency; final physical relation publication":"one physical local solve per authored occurrence"};
        if(contract.boundaryPreparation!=WorkspaceBoundary::None) {
            if(target.storageBinding!=StorageBinding::TransientWorkspace)
                throw std::runtime_error("LocalAffine boundary preparation requires its own transient workspace: step="+call.occurrence+" equation="+call.equation+" target="+target.id);
            result.requirements.push_back("nearest-interior extension of solved final coefficient, after internal evaluation; serial only, not MPI halo");
        }
        if(target.storageBinding!=StorageBinding::TransientWorkspace)
            result.physicalResources={fieldResource(call.occurrence,target.id,contract.valueReads,{},{},true)};
        if(target.storageBinding==StorageBinding::TransientWorkspace) {
            if(call.target.kind!=TargetKind::Physical || (target.role!=StateRole::Derived && target.role!=StateRole::Algebraic)
                || target.initializationRequired || target.boundaryRequired
                || target.outputEligible || target.restartEligible)
                throw std::runtime_error("Unsupported LocalAffine workspace lifecycle: "+target.id+"; no initial/boundary/output/restart port is implemented");
            result.requirements.push_back("solver-owned materialized workspace; no initial value; readable only after its HOW publication");
        }
        result.providerContract=std::move(contract);return result;
    }
namespace {
/// Block selection is a concrete steady numerical capability, independent of time recipes.
class LinearBlockMethod final : public IProvider {
public:
    std::string_view id() const override {return "LinearBlockCentral2Hypre";}
    std::string_view runtimeProvider() const override {return ScalarOps::Provider;}
    CompiledEquationCall compile(const ExecutableEquationSystem&,const EquationCall&,const NumericalBinding&) const override {
        throw std::runtime_error("LinearBlockCentral2Hypre requires a whole => equation/unknown set");
    }
    CompiledEquationCall compileBlock(const ExecutableEquationSystem& system,const ExecutionScope& scope,const NumericalBinding& binding) const override {
        if(scope.step.stageScoped)throw std::runtime_error("Time strategy conflict: implicit block inside explicit StageLoop is unsupported");
        if(scope.equations.size()!=scope.unknowns.size())throw std::runtime_error("Unsupported block layout: square scalar system required");
        if(!binding.inputs.empty() || !binding.operators.empty())throw std::runtime_error("Block layout comes from HOW; extra WHICH inputs/operators are unsupported by this fixed Central2 implementation");
        CompiledLinearBlock block;block.equations=scope.equations;block.targets=scope.unknowns;
        for(const auto& name:block.targets){const auto& s=system.state.at(name);if(s.readOnly() || s.producer || s.constantValue || s.components!=1 || s.shape!=ValueShape::Scalar || s.componentOffset!=0 || s.location!=VariableLocation::EulerianCell || (s.storageBinding!=StorageBinding::NamedDistributed && s.storageBinding!=StorageBinding::ProviderDistributed))throw std::runtime_error("Unsupported writable scalar block target layout: "+name);}
        for(const auto& parameter:binding.parameters){double v=parameter.second;if(!std::isfinite(v))throw std::runtime_error("Nonfinite block linear control");if(parameter.first=="relativeTolerance")block.linear.relativeTolerance=v;else if(parameter.first=="absoluteTolerance")block.linear.absoluteTolerance=v;else if(parameter.first=="maxIterations" || parameter.first=="krylovDimension"){if(v!=std::floor(v) || v<1 || v>1000000)throw std::runtime_error("Invalid block integer iteration control");if(parameter.first=="maxIterations")block.linear.maxIterations=v;else block.linear.krylovDimension=v;}else throw std::runtime_error("Unsupported block linear parameter: "+parameter.first);}
        FDM::validateLinearSolverConfig(block.linear,"LinearBlockCentral2Hypre");
        std::vector<std::string> symbols=block.targets;
        for(const auto& eq:block.equations){const auto& e=system.registry.at(eq);system.state.validateSymbols(e.lhs);system.state.validateSymbols(e.rhs);for(const auto* side:{&e.lhs,&e.rhs})for(const auto& name:inspectFormula(*side).symbols){const auto& state=system.state.at(name);if(!state.constantValue && std::find(symbols.begin(),symbols.end(),name)==symbols.end())symbols.push_back(name);}}
        for(const auto& name:symbols){const auto& s=system.state.at(name);if(s.components!=1 || s.shape!=ValueShape::Scalar || s.componentOffset!=0 || s.location!=VariableLocation::EulerianCell || (s.storageBinding!=StorageBinding::NamedDistributed && s.storageBinding!=StorageBinding::ProviderDistributed) || (s.producer && s.evaluation!=StateEvaluation::Lazy))throw std::runtime_error("Unsupported Current block STATE input layout: "+name);block.valueReads.push_back({name,StateVersion::Current});}
        for(const auto& read:scope.step.inputs)if(std::find(block.targets.begin(),block.targets.end(),read.symbol)!=block.targets.end() || read.kind!=TargetKind::Physical || read.version!=StateVersion::Current || read.halo || !read.snapshot.empty() || read.perIteration || std::none_of(block.valueReads.begin(),block.valueReads.end(),[&](const auto& r){return r.symbol==read.symbol;}))throw std::runtime_error("Unsupported block input version/publication: "+read.symbol);
        const auto depends=[&](const FormulaExpr& e){const auto names=inspectFormula(e).symbols;return std::any_of(block.targets.begin(),block.targets.end(),[&](const auto& n){return std::find(names.begin(),names.end(),n)!=names.end();});};
        const auto value=[&](const Equation& eq,const FormulaExpr& e){std::vector<std::string> used;return compileFormulaValue(eq,e,{},used,symbols,&system.state);};
        for(const auto& eqid:block.equations){
            const auto& eq=system.registry.at(eqid);std::string rowTarget;int diffusionCount=0;
            // Structural affine assembly relative to the entire unknown set; never sample at 0/1 or freeze nonlinear terms.
            const auto lower=[&](const auto& self,const FormulaExpr& e)->FormulaLinearKernel {
                using K=FormulaExpr::Kind;
                if(e.kind()==K::Operator){
                    if(e.operation()==Math::Op::Ddt)throw std::runtime_error("Time strategy conflict: steady block has no implicit time implementation; configured explicit recipe cannot advance =>");
                    if(e.operation()!=Math::Op::Diffusion || e.arguments().size()!=2 || e.arguments()[1].kind()!=K::Symbol)throw std::runtime_error("Unsupported block numerical operator");
                    const auto& input=e.arguments()[1].name();auto target=std::find(block.targets.begin(),block.targets.end(),input);if(target==block.targets.end() || depends(e.arguments()[0]))throw std::runtime_error("Unsupported block diffusion input/coefficient");
                    const auto& d=e.arguments()[0];double D;
                    if(d.kind()==K::Constant)D=d.constant();else if(d.kind()==K::Symbol && system.state.at(d.name()).constantValue)D=*system.state.at(d.name()).constantValue;else throw std::runtime_error("Block diffusion requires a declared constant coefficient");
                    if(!std::isfinite(D) || D<=0)throw std::runtime_error("Block diffusion coefficient must be positive and finite");
                    if(++diffusionCount!=1)throw std::runtime_error("Unsupported block row: exactly one principal diffusion required");rowTarget=input;
                    return [D,input](const auto& v,int cell,int){return Discretization::nodalCentral2Linear(D,input,v,cell);};
                }
                if(!depends(e)){auto f=value(eq,e);return [f](const auto& v,int cell,int c){FormulaLinearTerm row;row.constant=f(v,cell,c);return row;};}
                if(e.kind()==K::Symbol){auto name=e.name();return [name](const auto& v,int cell,int){FormulaLinearTerm row;row.add(v.dof(name,cell,0),1);return row;};}
                if(e.kind()==K::Negate){auto f=self(self,e.arguments()[0]);return [f](const auto& v,int cell,int c){auto row=f(v,cell,c);row.constant=-row.constant;for(auto& a:row.coefficients)a.value=-a.value;return row;};}
                if(e.kind()==K::Add || e.kind()==K::Subtract){auto l=self(self,e.arguments()[0]),r=self(self,e.arguments()[1]);double sign=e.kind()==K::Add?1:-1;return [l,r,sign](const auto& v,int cell,int c){auto row=l(v,cell,c),other=r(v,cell,c);row.constant+=sign*other.constant;for(const auto& a:other.coefficients)row.add(a.column,sign*a.value);return row;};}
                bool left=depends(e.arguments().at(0)),right=e.arguments().size()>1 && depends(e.arguments()[1]);
                if(e.kind()==K::Multiply || e.kind()==K::Divide){if((e.kind()==K::Multiply && left && right) || (e.kind()==K::Divide && right))throw std::runtime_error("Unsupported nonlinear block residual/target denominator");auto f=self(self,e.arguments()[left?0:1]);auto factor=value(eq,e.arguments()[left?1:0]);bool divide=e.kind()==K::Divide;return [f,factor,divide](const auto& v,int cell,int c){double q=factor(v,cell,c);if(!std::isfinite(q) || (divide && q==0))throw std::runtime_error("Invalid block coefficient/denominator");if(divide)q=1/q;auto row=f(v,cell,c);row.constant*=q;for(auto& a:row.coefficients)a.value*=q;return row;};}
                throw std::runtime_error("Unsupported nonlinear block residual");
            };
            try{auto residual=FormulaExpr::subtract(eq.lhs,eq.rhs);block.rows.push_back(lower(lower,residual));}catch(const std::exception& e){throw std::runtime_error("equation="+eqid+": "+e.what());}
            if(diffusionCount!=1 || std::find(block.rowTargets.begin(),block.rowTargets.end(),rowTarget)!=block.rowTargets.end())throw std::runtime_error("Block layout lacks unique principal row for every unknown");block.rowTargets.push_back(rowTarget);
        }
        CompiledEquationCall result;result.source=scope.step;result.target.symbol=scope.step.target.symbol;result.equationMethod=id();result.writes=block.targets;result.stateUses=block.valueReads;
        for(const auto& read:block.valueReads)result.reads.push_back(read.symbol);
        for(std::size_t i=0;i<block.equations.size();++i)result.calls.push_back({block.equations[i],block.rowTargets[i]});
        for(const auto& name:block.targets){result.stateEffects.push_back({name,true,false,true});result.physicalResources.push_back(fieldResource(instanceOperation(instanceOperation("block.data",scope.id),name),name,block.valueReads,block.targets));}
        result.fragment.kind=PlanNodeKind::Sequence;result.fragment.id=scope.id;result.fragment.name="steady scalar whole-system solve";result.fragment.equations=scope.equations;
        for(const auto& op:{"block.assemble","block.solve","block.publish"}){auto id=instanceOperation(op,scope.id);SolvePlanNode n;n.kind=std::string(op)=="block.assemble"?PlanNodeKind::Assemble:std::string(op)=="block.solve"?PlanNodeKind::BlockSolve:PlanNodeKind::Commit;n.id=n.name=n.operation=id;n.equationCalls=result.calls;result.fragment.children.push_back(n);result.operations.push_back({id,op,std::string(op)=="block.publish"?OperationStage::StepCommit:OperationStage::Prepare,{std::string(op)=="block.solve"?OperationCapability::SparseLinearSolve:OperationCapability::ScalarExplicit},{OriginKind::Generated,std::string(this->id())}});}
        result.requirements={"serial single patch; square scalar steady system","one principal constant Central2 diffusion per row; explicit HYPRE FlexGMRES/BoomerAMG","all residuals accepted before atomic multi-target publication; no physical TimeCommit"};result.providerContract=std::move(block);return result;
    }
};

}
void addScalarMethod(ProviderRegistry& registry) {
    static const ScalarMethod method,transport(ScalarTransportForm::Conservative),
        shared(ScalarTransportForm::SharedConservative),stage(ScalarTransportForm::SharedPrimitive),
        variable(ScalarTransportForm::None,true);
    for(const auto* provider:{&method,&transport,&shared,&stage,&variable})registry.add(*provider);
    static const LinearBlockMethod block;registry.add(block);
}
}
