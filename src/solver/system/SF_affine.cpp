#include "SF_affine.h"
#include "SF_formulaCompiler.h"
#include "SF_stateRealization.h"
#include "core/system/SF_operationIds.h"

namespace SF::System {
/// @brief 精确分解局部仿射式，系数来自表达式结构。
CompiledEquationCall compileLocalAffine(const ExecutableEquationSystem& system,const EquationCall& call,const NumericalBinding& binding) {
    if(!binding.inputs.empty() || !binding.parameters.empty() || !binding.operators.empty())throw std::runtime_error("LocalAffine accepts STATE read versions, not hidden method/operator parameters");
    const auto& target=system.state.at(call.target.symbol);
    if(call.target.kind==TargetKind::Workspace || (call.stageScoped && call.target.kind!=TargetKind::Physical) || target.readOnly() || target.shape!=ValueShape::Scalar || target.components!=1 || target.componentOffset!=0 || target.location!=VariableLocation::EulerianCell || (target.storageBinding!=StorageBinding::NamedDistributed && target.storageBinding!=StorageBinding::ProviderDistributed && target.storageBinding!=StorageBinding::TransientWorkspace))throw std::runtime_error("LocalAffine requires a writable scalar cell target with an existing or solver workspace binding");
    const auto& equation=system.registry.at(call.equation);const auto residual=FormulaExpr::subtract(equation.lhs,equation.rhs);const auto shape=inspectFormula(residual,call.target.symbol);
    if(std::find(shape.symbols.begin(),shape.symbols.end(),call.target.symbol)==shape.symbols.end())throw std::runtime_error("HOW target is absent from equation: "+call.target.symbol);
    if(shape.targetDegree!=1)throw std::runtime_error("Unsupported LocalAffine: nonlinear or target-dependent denominator; select an implemented nonlinear method with branch/convergence control");
    for(auto operation:shape.operations)if(operation!=Math::Op::Add && operation!=Math::Op::Subtract && operation!=Math::Op::Multiply && operation!=Math::Op::Divide && operation!=Math::Op::Negate && operation!=Math::Op::Pow && operation!=Math::Op::Component && operation!=Math::Op::Vector)throw std::runtime_error("LocalAffine supports local algebra only; no differential/operator solve is implemented");
    CompiledAffineEquation contract;contract.boundaryPreparation=target.workspaceBoundary;
    std::vector<StateUse> valueReads;
    std::vector<std::string> slots;
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
        valueReads.push_back(use);slots.push_back(name);
    }
    for(const auto& read:call.inputs)
        if(std::count_if(call.inputs.begin(),call.inputs.end(),[&](const auto& r){return r.symbol==read.symbol;})!=1
            || std::none_of(valueReads.begin(),valueReads.end(),[&](const auto& r){return r.symbol==read.symbol;}))
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
    CompiledEquationCall result;result.source=call;result.target.symbol=target.id;result.target.kind=call.target.kind;result.stageLocal=call.stageScoped;result.equationMethod="LocalAffine";result.backendProvider=ScalarOps::Provider;result.calls={{call.equation,target.id}};result.writes={target.id};result.stateUses=valueReads;for(const auto& input:valueReads)result.reads.push_back(input.symbol);result.stateEffects={{target.id,true,false,true,true}};
    if(call.target.kind!=TargetKind::Physical) {
        result.target.workspace="state."+target.id+(call.target.kind==TargetKind::Working?".Working":".Correction");
        result.stateEffects.clear();result.workspaceProvides={result.target.workspace};
        result.stepWorkspaces=result.workspaceProvides;
    }
    for(const auto& input:valueReads)if(input.kind!=TargetKind::Physical)
        result.workspaceRequires.push_back("state."+input.symbol+(input.kind==TargetKind::Working?".Working":".Correction"));
    result.residualCapabilities={ResidualKind::StateChange,ResidualKind::DiscreteEquation};
    result.backendOperation=instanceOperation("affine.solve",call.occurrence);result.operations.push_back({result.backendOperation,"solve local affine residual for "+target.id,OperationStage::Prepare,{OperationCapability::ScalarExplicit},{OriginKind::Generated,"LocalAffine"}});
    if(result.stageLocal){result.stageWorkspace="affine."+call.occurrence+".stage";result.numericalWorkspaces={result.stageWorkspace};result.temporalPublish=instanceOperation("affine.publish",call.occurrence);result.operations.push_back({result.temporalPublish,"publish algebraic relation from committed physical inputs",OperationStage::StepCommit,{OperationCapability::ScalarExplicit},{OriginKind::Generated,"LocalAffine"}});}
    result.requirements={"serial single patch","local scalar affine residual; no root branch selection",result.stageLocal?"explicit HOW Stage frequency; final physical relation publication":"one local solve per authored invocation; no clock ownership"};
    if(contract.boundaryPreparation!=WorkspaceBoundary::None) {
        if(target.storageBinding!=StorageBinding::TransientWorkspace)
            throw std::runtime_error("LocalAffine boundary preparation requires its own transient workspace: step="+call.occurrence+" equation="+call.equation+" target="+target.id);
        result.requirements.push_back("nearest-interior extension of solved final coefficient, after internal evaluation; serial only, not MPI halo");
    }
    if(target.storageBinding!=StorageBinding::TransientWorkspace)
        result.physicalResources={equationFieldResource(call.occurrence,target.id,valueReads,{},{},true)};
    if(target.storageBinding==StorageBinding::TransientWorkspace) {
        if(call.target.kind!=TargetKind::Physical || (target.role!=StateRole::Derived && target.role!=StateRole::Algebraic)
            || target.initializationRequired || target.boundaryRequired
            || target.outputEligible || target.restartEligible)
            throw std::runtime_error("Unsupported LocalAffine workspace lifecycle: "+target.id+"; no initial/boundary/output/restart port is implemented");
        result.requirements.push_back("solver-owned materialized workspace; no initial value; readable only after its HOW publication");
    }
    result.providerContract=std::move(contract);return result;
}
} // namespace SF::System
