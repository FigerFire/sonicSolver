#include "SF_cell.h"
#include "SF_formulaCompiler.h"
#include "SF_stateRealization.h"
#include "SF_numericalSystem.h"
#include "SF_termProviderCatalog.h"
#include "core/system/SF_operationIds.h"
#include <cmath>

namespace SF::System {
namespace {
/// @brief 物化系数和 RHS lowering 共用的应力散度形状。
bool isStressDivergence(const FormulaExpr& expression) {
    return expression.operation()==Math::Op::Div && expression.arguments().size()==1
        && expression.arguments()[0].kind()==FormulaExpr::Kind::Multiply
        && expression.arguments()[0].arguments()[1].operation()==Math::Op::DevSymmGrad;
}
}
CompiledEquationCall compileCellEquation(const ExecutableEquationSystem& system,const EquationCall& call,const NumericalBinding& binding) {
    const auto& target=system.state.at(call.target.symbol);
    if((call.target.kind!=TargetKind::Physical && call.target.kind!=TargetKind::Working) || target.components!=3
        || target.shape!=ValueShape::Vector || target.readOnly() || target.location!=VariableLocation::EulerianCell
        || (target.storageBinding!=StorageBinding::PackedDistributed && target.storageBinding!=StorageBinding::NamedDistributed))
        throw std::runtime_error("CellExplicit requires writable three-component cell storage");
    if(!binding.inputs.empty() || !binding.parameters.empty() || !binding.response.empty())
        throw std::runtime_error("CellExplicit does not consume extra WHICH inputs, parameters or response references");
    const auto& equation=system.registry.at(call.equation);
    explicitFormulaRhs(equation,target.id); // The existing strict ddt(target) contract remains authoritative.
    CompiledEquationCall result;result.source=call;result.target.symbol=target.id;result.target.kind=call.target.kind;
    const auto work="cell."+call.occurrence+".working";
    if(call.target.kind==TargetKind::Working){result.target.workspace=work;result.target.viewOwner=StateViewOwner::NumericalProvider;}
    result.equationMethod="CellExplicit";result.backendProvider="equation.cell-explicit";result.spatialTerms=true;
    result.temporalCapabilities={FDM::TimeRecipeId::ForwardEuler,FDM::TimeRecipeId::PartitionedFixedEuler};result.calls={{equation.id,target.id}};result.writes={target.id};
    CompiledStateView view;view.symbol=target.id;view.kind=StateViewKind::Working;view.storage=work;
    view.components=3;view.owner=StateViewOwner::NumericalProvider;view.haloRequired=true;result.viewOutputs={view};
    result.workspaceProvides={"cell.response:"+call.occurrence};result.stepWorkspaces=result.workspaceProvides;
    result.workspaceInvalidates=result.workspaceProvides;
    result.stepInputs={{OpIds::CellPrepare,"equation.cell-explicit",true}};
    result.physicalPublication={OpIds::CellAccept,"equation.cell-explicit"};
    result.providerContract=CompiledCellEquation{};
    auto port=equationFieldResource("cell.field."+call.occurrence,target.id,{},{});
    const auto spatial=[&](const auto& self,const FormulaExpr& expression)->bool {
        if(expression.operation()==Math::Op::Div || expression.operation()==Math::Op::Diffusion || expression.operation()==Math::Op::Grad)return true;
        return std::any_of(expression.arguments().begin(),expression.arguments().end(),[&](const auto& child){return self(self,child);});
    };
    port.fields->boundaryRequired=spatial(spatial,equation.rhs);port.fields->localOnly=!port.fields->boundaryRequired;port.fields->cellVolumeBoundary=true;
    if(port.fields->boundaryRequired)port.fields->boundaries.push_back(target.id);
    result.physicalResources.push_back(std::move(port));
    if(call.target.kind==TargetKind::Physical)result.stateEffects={{target.id,true,true}};
    for(const auto& id:inspectFormula(equation.rhs).symbols)if(system.state.contains(id) && !system.state.at(id).constantValue) {
        const auto selected=std::find_if(call.inputs.begin(),call.inputs.end(),[&](const auto& r){return r.symbol==id;});
        StateUse use{id,call.stageScoped?StateVersion::Stage:StateVersion::Current};if(selected!=call.inputs.end())use=*selected;
        const bool frozen=use.version==StateVersion::Frozen && use.perIteration && !call.stageScoped && id!=target.id;
        const bool working=use.kind==TargetKind::Working && use.version==StateVersion::Current && id!=target.id;
        if((use.version!=StateVersion::Current && !frozen && !(call.stageScoped && (use.version==StateVersion::Stage || use.version==StateVersion::OldTime))) || (use.kind!=TargetKind::Physical && !working) || use.halo || !use.snapshot.empty() || (use.perIteration && !frozen))
            throw std::runtime_error("CellExplicit input="+id+" requires implemented Current (Step) or Stage/OldTime (recipe) backing and unqualified access");
        result.stateUses.push_back(use);result.reads.push_back(id);
    }
    if(std::none_of(result.stateUses.begin(),result.stateUses.end(),[&](const auto& r){return r.symbol==target.id;}))result.stateUses.push_back({target.id,call.stageScoped?StateVersion::Stage:StateVersion::Current});
    for(const auto& use:call.inputs)if(std::count_if(call.inputs.begin(),call.inputs.end(),[&](const auto& other){return other.symbol==use.symbol;})!=1
        || std::none_of(result.stateUses.begin(),result.stateUses.end(),[&](const auto& input){return input.symbol==use.symbol;}))
        throw std::runtime_error("Duplicate/unconsumed cell STATE input: "+call.occurrence+" / "+use.symbol);
    result.physicalResources.front().fields->states.clear();
    result.physicalResources.front().fields->states.push_back(target.id);
    for(const auto& use:result.stateUses)result.physicalResources.front().fields->states.push_back(use.symbol);
    const auto inputs=[&](const auto& self,const FormulaExpr& expression)->void {
        if(expression.operation()==Math::Op::DevSymmGrad && expression.arguments().size()==1 && expression.arguments()[0].kind()==FormulaExpr::Kind::Symbol && expression.arguments()[0].name()!=target.id)
            throw std::runtime_error("Cell stress currently acts on this call target vector: "+call.occurrence);
        if(isStressDivergence(expression)) {
            const auto& coefficient=expression.arguments()[0].arguments()[0];
            if(coefficient.kind()==FormulaExpr::Kind::Symbol) {
                const auto& symbol=system.state.at(coefficient.name());
                if(!symbol.constantValue) {
                    if(symbol.components!=1 || symbol.componentOffset || symbol.location!=VariableLocation::EulerianCell || symbol.producer)
                        throw std::runtime_error("Cell stress requires an unsliced materialized cell coefficient: "+symbol.id);
                    const auto use=std::find_if(result.stateUses.begin(),result.stateUses.end(),[&](const auto& read){return read.symbol==symbol.id;});
                    if(use==result.stateUses.end() || use->version==StateVersion::OldTime)
                        throw std::runtime_error("Cell stress coefficient needs an implemented spatial view: "+symbol.id);
                    if(symbol.storageBinding==StorageBinding::TransientWorkspace) {
                        if(symbol.workspaceBoundary!=WorkspaceBoundary::NearestInterior)throw std::runtime_error("Cell stress coefficient needs an explicit spatial publication: "+symbol.id);
                    } else if(use->version==StateVersion::Current && use->kind==TargetKind::Physical) {
                        auto input=equationFieldResource("cell.input."+call.occurrence+"."+symbol.id,symbol.id,{});
                        input.fields->boundaryRequired=true;input.fields->cellVolumeBoundary=true;input.fields->boundaries.push_back(symbol.id);
                        if(std::none_of(result.physicalResources.begin(),result.physicalResources.end(),[&](const auto& port){return port.id==input.id;}))result.physicalResources.push_back(std::move(input));
                    }
                }
            }
        }
        if(expression.operation()==Math::Op::Grad && expression.arguments().size()==1 && expression.arguments()[0].kind()==FormulaExpr::Kind::Symbol) {
            const auto& id=expression.arguments()[0].name();
            const auto use=std::find_if(result.stateUses.begin(),result.stateUses.end(),[&](const auto& value){return value.symbol==id;});
            if(call.stageScoped && (use==result.stateUses.end() || use->version!=StateVersion::Stage))throw std::runtime_error("Cell Stage gradient requires a real Stage boundary port: "+call.occurrence+" / "+id);
            auto input=equationFieldResource("cell.input."+call.occurrence+"."+id,id,{});
            input.fields->boundaryRequired=true;input.fields->cellVolumeBoundary=true;input.fields->boundaries.push_back(id);
            if(std::none_of(result.physicalResources.begin(),result.physicalResources.end(),[&](const auto& port){return port.id==input.id;}))result.physicalResources.push_back(std::move(input));
        }
        for(const auto& arg:expression.arguments())self(self,arg);
    };inputs(inputs,equation.rhs);
    for(const auto* op:{OpIds::CellAssemble,OpIds::CellAdvance}) {
        SolvePlanNode leaf;leaf.kind=std::string_view(op)==OpIds::CellAssemble?PlanNodeKind::Assemble:PlanNodeKind::Update;
        leaf.id=leaf.name=leaf.operation=instanceOperation(op,call.occurrence);leaf.provider=result.backendProvider;leaf.equationCalls=result.calls;result.fragment.children.push_back(leaf);
        result.operations.push_back({leaf.operation,op,OperationStage::Prepare,{OperationCapability::ScalarExplicit},{OriginKind::Generated,result.equationMethod}});
    }
    for(const auto* op:{OpIds::CellPrepare,OpIds::CellAccept})result.operations.push_back({op,op,OperationStage::Prepare,{OperationCapability::ScalarExplicit},{OriginKind::Generated,result.equationMethod}});
    if(call.stageScoped) {
        if(call.target.kind!=TargetKind::Physical)throw std::runtime_error("Cell recipe participant requires Physical acceptance target: "+call.occurrence);
        result.temporalResidual=result.synchronousStages=true;result.temporalStepSize=TemporalOps::Dt;
        result.temporalCapabilities={FDM::TimeRecipeId::ForwardEuler,FDM::TimeRecipeId::SSPRK3,FDM::TimeRecipeId::ClassicalRK4};
        result.stepInputs.clear();result.physicalPublication={};result.workspaceProvides.clear();result.stepWorkspaces.clear();result.workspaceInvalidates.clear();
        result.fragment={};result.operations.clear();result.viewOutputs.clear();
        result.oldTimeWorkspace="cell."+call.occurrence+".old";result.stageWorkspace="cell."+call.occurrence+".stage";result.residualWorkspace="cell."+call.occurrence+".rhs";
        result.temporalSnapshot=instanceOperation("cell.snapshot",call.occurrence);
        result.temporalStagePrepare=instanceOperation("cell.prepare",call.occurrence);
        result.temporalRhs=instanceOperation("cell.rhs",call.occurrence);
        result.temporalAdvance=instanceOperation("cell.advance",call.occurrence);
        result.temporalPublish=instanceOperation("cell.publish",call.occurrence);
        result.stateEffects={{result.residualWorkspace,true,false,false}};
        for(const auto& op:{result.temporalSnapshot,result.temporalStagePrepare,result.temporalRhs,result.temporalAdvance,result.temporalPublish})
            result.operations.push_back({op,op,OperationStage::Prepare,{OperationCapability::ScalarExplicit},{OriginKind::Generated,result.equationMethod}});
        for(const auto& use:result.stateUses)if(use.symbol==target.id && use.version!=StateVersion::Stage)
            throw std::runtime_error("Cell integrated stencil requires its real Stage input: "+call.occurrence);
    }
    result.requirements={"vector additive spatial and point contributions, strict ddt(target)","per-occurrence workspace and explicit Euler response; common physical interval, optional explicitly bound fixed base","no physical model/equation identity or automatic algorithm selection"};
    return result;
}

void lowerCellEquation(const ExecutableEquationSystem& system,CompiledEquationCall& call,
        const FDM::NumericalRecipeSet& recipes,const TermProviderCatalog& catalog,
        const std::vector<FormulaOperatorBinding>& bindings,CompiledNumericalSystem& result) {
    const auto& equation=system.registry.at(call.source.equation);
    const auto firstTerm=result.operators.size();
    const auto fail=[&](const std::string& occurrence,const std::string& why) {
        throw std::runtime_error("Cell equation="+equation.id+" step="+call.source.occurrence+" target="+call.target.symbol+" occurrence="+occurrence+": "+why);
    };
    CompiledCellEquation contract;
    const auto& lhs=equation.lhs;
    if(lhs.kind()!=FormulaExpr::Kind::Operator || lhs.operation()!=Math::Op::Ddt || lhs.arguments().size()!=1
        || lhs.arguments()[0].kind()!=FormulaExpr::Kind::Symbol || lhs.arguments()[0].name()!=call.target.symbol)
        fail("lhs","this implementation requires ddt(target) alone on LHS; no implicit mass normalization");
    const auto timeOccurrence=lhs.occurrence().empty()?"lhs":lhs.occurrence();
    if(selectFormulaBinding(equation,lhs,timeOccurrence,bindings))fail(timeOccurrence,"time uses the unique recipe, not a spatial override");
    const auto constant=[&](const FormulaExpr& e)->double {
        if(e.kind()==FormulaExpr::Kind::Constant)return e.constant();
        if(e.kind()==FormulaExpr::Kind::Symbol && system.state.at(e.name()).constantValue)return *system.state.at(e.name()).constantValue;
        fail(e.occurrence(),"coefficient needs a declared constant; no variable-field stencil is implemented here");return 0;
    };
    const auto visit=[&](const auto& self,const FormulaExpr& e,double sign,const std::string& path)->void {
        using K=FormulaExpr::Kind;const auto occurrence=e.occurrence().empty()?path:e.occurrence();
        if(!std::isfinite(sign))fail(occurrence,"nonfinite contribution coefficient");
        if(e.kind()==K::Add || e.kind()==K::Subtract) {
            self(self,e.arguments()[0],sign,path+"/0");self(self,e.arguments()[1],e.kind()==K::Add?sign:-sign,path+"/1");return;
        }
        if(e.kind()==K::Negate){self(self,e.arguments()[0],-sign,path+"/0");return;}
        if(e.kind()==K::Constant && e.constant()==0)return;
        if(e.kind()==K::Multiply) {
            const auto& a=e.arguments();
            if(a[0].kind()==K::Constant || (a[0].kind()==K::Symbol && system.state.at(a[0].name()).constantValue)) {self(self,a[1],sign*constant(a[0]),path+"/1");return;}
            fail(occurrence,"scaled contributions require a left constant coefficient");
        }
        const auto& stress=e.kind()==K::Divide?e.arguments()[0]:e;
        if(isStressDivergence(stress)) {
            const auto& product=stress.arguments()[0].arguments();
            if(product[1].arguments().size()!=1 || product[1].arguments()[0].kind()!=K::Symbol || product[1].arguments()[0].name()!=call.target.symbol)
                fail(occurrence,"deviatoric symmetric gradient must act on this vector target");
            const auto& coefficient=product[0];
            if(coefficient.kind()!=K::Constant && coefficient.kind()!=K::Symbol)fail(occurrence,"stress coefficient needs a constant or a published scalar; materialize coefficient relations explicitly");
            const double divisor=e.kind()==K::Divide?constant(e.arguments()[1]):1;
            if(!std::isfinite(divisor) || divisor<=0)fail(occurrence,"stress divisor needs a positive declared constant density");
            const auto address=stress.occurrence().empty()?(e.kind()==K::Divide?path+"/0":path):stress.occurrence();
            const auto* choice=selectFormulaBinding(equation,stress,address,bindings);
            if(!choice || choice->provider!="stress.central2.arithmetic-faces")fail(address,"select stress.central2.arithmetic-faces explicitly; no implicit face coefficient interpolation");
            const auto selected=catalog.resolve({equation,stress,address,call.target.symbol,call.equationMethod,nullptr,&system.state,&call.stateUses},bindings);
            if(selected.status!=BindingStatus::Resolved || !selected.primitiveSpatial)fail(address,"selected stress implementation has no bound kernel");
            contract.terms.push_back(result.operators.size());
            result.operators.emplace_back(equation.id,equation.id,call.target.symbol,call.equationMethod,address,Math::Op::Div,call.target.symbol,
                coefficient.kind()==K::Symbol?coefficient.name():"constant",CompiledSpatialBinding::Side::Right,std::nullopt,selected.id,selected.owner,true,PrimitiveMomentumSource{},ConservativeSourceKernel{},selected.primitiveSpatial);
            auto& bound=result.operators.back();bound.coefficient=sign/divisor;
            if(coefficient.kind()==K::Symbol && !system.state.at(coefficient.name()).constantValue) {
                bound.coefficientState=coefficient.name();
            }
            else {bound.constantCoefficient=constant(coefficient);if(bound.constantCoefficient<0)fail(address,"negative stress coefficient");contract.diffusivity+=4./3.*std::abs(bound.coefficient)*bound.constantCoefficient;}
            result.requiredHaloWidth=std::max(result.requiredHaloWidth,1);return;
        }
        const bool divided=e.kind()==K::Divide && e.arguments()[0].operation()==Math::Op::Grad;
        if(divided || e.operation()==Math::Op::Grad) {
            const auto& grad=divided?e.arguments()[0]:e;
            if(grad.arguments().size()!=1 || grad.arguments()[0].kind()!=K::Symbol)fail(occurrence,"gradient needs a registered scalar input");
            const auto& input=system.state.at(grad.arguments()[0].name());
            if(input.components!=1 || input.location!=VariableLocation::EulerianCell)fail(occurrence,"Central2 gradient requires a scalar cell input, not a face/constraint layout");
            double divisor=divided?constant(e.arguments()[1]):1.0;
            if(!std::isfinite(divisor) || divisor<=0)fail(occurrence,"gradient divisor must be positive and finite");
            const auto address=grad.occurrence().empty()?(divided?path+"/0":path):grad.occurrence();
            const auto selected=catalog.resolve({equation,grad,address,call.target.symbol,call.equationMethod,nullptr,&system.state,&call.stateUses},bindings);
            if(selected.status!=BindingStatus::Resolved || !selected.gradient)fail(address,selected.reason.empty()?"selected gradient has no actual kernel":selected.reason);
            contract.gradients.push_back({grad.arguments()[0].name(),sign,divisor,selected.gradient});
            if(divided && sign==-1 && contract.gradients.size()==1) {
                contract.potential=grad.arguments()[0].name();contract.capacity=e.arguments()[1].kind()==K::Symbol?e.arguments()[1].name():std::string{};contract.capacityValue=divisor;
            } else {contract.potential.clear();contract.capacity.clear();}
            result.operators.emplace_back(equation.id,equation.id,call.target.symbol,call.equationMethod,address,Math::Op::Grad,grad.arguments()[0].name(),divided?"constant divisor":"",CompiledSpatialBinding::Side::Right,std::nullopt,selected.id,selected.owner,true);
            result.operators.back().gradient=selected.gradient;result.operators.back().coefficient=sign/divisor;
            result.requiredHaloWidth=std::max(result.requiredHaloWidth,1);return;
        }
        if(e.kind()!=K::Operator)fail(occurrence,"term has no implemented spatial/source contribution");
        const FDM::TermRecipe* recipe=nullptr;std::string primary=call.target.symbol,secondary;
        if(e.operation()==Math::Op::Div) {
            if(e.arguments().size()!=1 || e.arguments()[0].operation()!=Math::Op::Outer)fail(occurrence,"convection requires div(outer(target,target))");
            const auto& a=e.arguments()[0].arguments();
            if(a.size()!=2 || a[0].kind()!=K::Symbol || a[1].kind()!=K::Symbol || a[0].name()!=call.target.symbol || a[1].name()!=call.target.symbol)
                fail(occurrence,"selected self-advection stencil needs the same registered vector operands");
            if(!recipes.convection || recipes.convection->id()!=FDM::TermRecipeId::PrimitiveUpwind1)fail(occurrence,"PrimitiveUpwind1 stencil required; other convection kernels not implemented");
            recipe=&*recipes.convection;
            contract.advectionWeight+=std::abs(sign);
        } else if(e.operation()==Math::Op::Diffusion) {
            if(e.arguments().size()!=2 || e.arguments()[1].kind()!=K::Symbol || e.arguments()[1].name()!=call.target.symbol)fail(occurrence,"diffusion needs this target vector");
            sign*=constant(e.arguments()[0]);contract.diffusivity+=std::abs(sign);
            if(!std::isfinite(sign))fail(occurrence,"nonfinite diffusion coefficient");
            if(!recipes.diffusion || recipes.diffusion->id()!=FDM::TermRecipeId::Central2Explicit)fail(occurrence,"Central2Explicit stencil required");
            recipe=&*recipes.diffusion;secondary=e.arguments()[0].kind()==K::Symbol?e.arguments()[0].name():"constant";
        } else if(e.operation()!=Math::Op::Source)fail(occurrence,"unimplemented operator "+std::string(Math::name(e.operation())));
        const auto selected=catalog.resolve({equation,e,occurrence,call.target.symbol,call.equationMethod,recipe,&system.state,&call.stateUses},bindings);
        if(selected.status!=BindingStatus::Resolved)fail(occurrence,selected.reason);
        if(recipe && !selected.primitiveSpatial)fail(occurrence,"selected stencil has no vector kernel");
        if(!recipe && !selected.primitiveSource && !selected.fieldSource.bind)fail(occurrence,"source has no actual vector reader/kernel");
        contract.terms.push_back(result.operators.size());
        result.operators.emplace_back(equation.id,equation.id,call.target.symbol,call.equationMethod,occurrence,e.operation(),primary,secondary,CompiledSpatialBinding::Side::Right,selected.recipe,selected.id,selected.owner,selected.compiledDataAvailable,selected.primitiveSource,ConservativeSourceKernel{},selected.primitiveSpatial);
        auto& bound=result.operators.back();bound.coefficient=sign;bound.fieldSource=selected.fieldSource;
        if(selected.recipe)result.recipeBindings.push_back({*selected.recipe,RecipeConsumerKind::EquationTerm,equation.id+"@"+occurrence});
        if(recipe)result.requiredHaloWidth=std::max(result.requiredHaloWidth,recipe->haloWidth());
    };
    visit(visit,equation.rhs,1,"rhs");
    for(const auto& binding:bindings)if(binding.formula==equation.id || (binding.formula.empty() && binding.operation)) {
        if(std::none_of(result.operators.begin()+firstTerm,result.operators.end(),[&](const auto& term){return term.formulaId==equation.id
            && (binding.operation?*binding.operation==term.operation:binding.occurrence==term.occurrence);}))
            fail(binding.occurrence,"unconsumed operator choice "+binding.provider);
    }
    if(contract.gradients.size()!=1){contract.potential.clear();contract.capacity.clear();}
    call.providerContract=std::move(contract);
}

} // namespace SF::System
