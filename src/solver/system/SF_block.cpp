#include "SF_block.h"
#include "SF_formulaCompiler.h"
#include "SF_stateRealization.h"
#include "core/system/SF_operationIds.h"
#include "solver/discretization/diffusion/SF_formulaCentral2.h"
#include <cmath>

namespace SF::System {
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
        for(const auto& name:block.targets){result.stateEffects.push_back({name,true,false,true});result.physicalResources.push_back(equationFieldResource(instanceOperation(instanceOperation("block.data",scope.id),name),name,block.valueReads,block.targets));}
        result.fragment.kind=PlanNodeKind::Sequence;result.fragment.id=scope.id;result.fragment.name="steady scalar whole-system solve";result.fragment.equations=scope.equations;
        for(const auto& op:{"block.assemble","block.solve","block.publish"}){auto id=instanceOperation(op,scope.id);SolvePlanNode n;n.kind=std::string(op)=="block.assemble"?PlanNodeKind::Assemble:std::string(op)=="block.solve"?PlanNodeKind::BlockSolve:PlanNodeKind::Commit;n.id=n.name=n.operation=id;n.equationCalls=result.calls;result.fragment.children.push_back(n);result.operations.push_back({id,op,std::string(op)=="block.publish"?OperationStage::StepCommit:OperationStage::Prepare,{std::string(op)=="block.solve"?OperationCapability::SparseLinearSolve:OperationCapability::ScalarExplicit},{OriginKind::Generated,std::string(this->id())}});}
        result.requirements={"serial single patch; square scalar steady system","one principal constant Central2 diffusion per row; explicit HYPRE FlexGMRES/BoomerAMG","all residuals accepted before atomic multi-target publication; no physical TimeCommit"};result.providerContract=std::move(block);return result;
    }
};

} // namespace

const IProvider& linearBlockMethod() {
    static const LinearBlockMethod method;
    return method;
}
} // namespace SF::System
