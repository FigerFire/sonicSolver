/// @file SF_numericalCompiler.cpp
/// @brief Bind selected Equation AST occurrences to spatial numerical providers.

#include "SF_numericalCompiler.h"
#include "SF_config.h"
#include "SF_numericalSelection.h"

#include <algorithm>
#include <stdexcept>

namespace SF::System::NumericalCompiler {
namespace {

void appendUnique(std::vector<std::string>& values,std::string value) {
    if (std::find(values.begin(),values.end(),value)==values.end())
        values.push_back(std::move(value));
}

} // namespace

void resolveConservativeInputs(const FDM::SolverConfig& config,CompiledNumericalSystem& system) {
    config.validateThermophysicalProjection();
    const auto* thermo=config.thermophysical.get();
    system.thermophysical=config.thermophysical;
    system.conservative=ConservativeExecutionInputs{config.numerics.ibmBoundary,config.numerics.interfaceFlux,config.numerics.ilwOrder,
        thermo?thermo->dynamicViscosity:config.numerics.dynamicViscosity,
        thermo?thermo->prandtl:config.numerics.prandtl,
        thermo?thermo->gamma:config.numerics.idealGasGamma,
        thermo?thermo->gasConstant:config.numerics.idealGasConstant,
        thermo?thermo->source:"explicit low-level numerics parameters"};
}
const ConservativeExecutionInputs& requireConservativeInputs(const CompiledNumericalSystem& system) {
    if(!system.conservative)throw std::runtime_error("Conservative execution inputs were not resolved at composition");
    return *system.conservative;
}

CompiledEquationCall compileCall(
        const ExecutableEquationSystem& system,const StateRegistry& state,
        const ExecutionScope& site,const std::vector<NumericalBinding>& bindings,
        const ProviderRegistry& registry,const CompiledTimeRecipe* time) {
    const auto& step=site.step;
    const bool block=site.kind==ExecutionKind::BlockSolve;
    auto effective=resolveNumerics(step,bindings,block,&system);
    const auto& selected=effective.value;
    const auto& mathematical=system.registry.at(step.equation);
    const IProvider* legacy=nullptr;
    CompiledEquationCall compiled;
    try {
        state.validateSymbols(mathematical.lhs);state.validateSymbols(mathematical.rhs);
        auto direct=block?std::nullopt:compileDirectEquation(system,step,selected);
        if(direct)compiled=std::move(*direct);
        else {
            legacy=&registry.at(selected.method);
            compiled=block?legacy->compileBlock(system,site,selected):legacy->compile(system,step,selected);
            compiled.lowerSpatial=legacy->spatialLowering();
        }
    }
    catch(const std::exception& error) {
        throw std::runtime_error(mathematical.origin.source+": occurrence="+step.occurrence+" equation="+step.equation+" target="+step.target.symbol+" provider="+selected.method+": "+error.what());
    }
    for(const auto& read:step.inputs)
        if(std::none_of(compiled.stateUses.begin(),compiled.stateUses.end(),[&](const auto& use){return use.symbol==read.symbol && use.version==read.version
                && use.halo==read.halo && use.snapshot==read.snapshot && use.perIteration==read.perIteration && use.kind==read.kind;}))
            throw std::runtime_error("Unconsumed STATE input version: "+read.symbol+" at "+step.occurrence);
    if(!compiled.lowerSpatial)for(const auto& op:selected.operators)
        if(std::find(compiled.operatorBindings.begin(),compiled.operatorBindings.end(),op.provider)==compiled.operatorBindings.end())
            throw std::runtime_error("Unconsumed numerical operator choice: occurrence="+step.occurrence+" equation="+step.equation+" operator="+op.occurrence+" provider="+op.provider);

    if(selected.time && (!compiled.temporalResidual || !time || *selected.time!=time->id()))
        throw std::runtime_error("Time recipe conflict/unused ddt choice: occurrence="+step.occurrence+" equation="+step.equation+" provider="+selected.method+"; one resolved time recipe is required");
    compiled.legacyProvider=legacy!=nullptr;
    if(legacy) {
        compiled.backendProvider=legacy->runtimeProvider();
        const auto ownFragment=[&](const auto& self,SolvePlanNode& node)->void {
            if(!node.operation.empty())node.provider=legacy->operationProvider(node.operation);
            for(auto& child:node.children)self(self,child);
        };
        ownFragment(ownFragment,compiled.fragment);
    }
    if (compiled.source.equation!=step.equation || compiled.source.occurrence!=step.occurrence
        || compiled.source.target.symbol!=step.target.symbol || compiled.source.target.kind!=step.target.kind
        || compiled.target.symbol!=step.target.symbol || compiled.target.kind!=step.target.kind)
        throw std::runtime_error("Provider changed frozen equation occurrence/target semantics: "+step.occurrence);
    if (time && !compiled.temporalCapabilities.empty()
        && std::find(compiled.temporalCapabilities.begin(),compiled.temporalCapabilities.end(),time->id())
            ==compiled.temporalCapabilities.end())
        throw std::runtime_error("Unsupported: provider "+selected.method+" cannot compile occurrence "
            +step.occurrence+" with selected temporal method "+FDM::toString(time->id()));
    if (!compiled.fragment.children.empty() || !compiled.fragment.operation.empty()) {
        if (!compiled.backendOperation.empty())
            throw std::runtime_error("EquationMethod has both fragment and backend operation authorities.");
        const auto validate=[&](const auto& self,const SolvePlanNode& node)
            -> void {
            if (node.children.empty() && node.operation.empty())
                throw std::runtime_error("EquationMethod fragment has a leaf without OpId.");
            for (const auto& child:node.children) self(self,child);
        };
        validate(validate,compiled.fragment);
    } else if (compiled.backendOperation.empty()) {
        throw std::runtime_error("EquationMethod supplied neither a fragment nor backend OpId.");
    }
    if (compiled.stateEffects.empty() && compiled.temporalResidual
        && compiled.publishesStageToPhysicalTarget && compiled.target.kind==TargetKind::Physical) {
        for (const auto& written:compiled.writes) if (state.contains(written))
            compiled.stateEffects.push_back({written,false,false});
    }
    compiled.numericalOperators=std::move(effective.value.operators);
    compiled.numericalOrigins=std::move(effective.origins);
    return compiled;
}

CompiledNumericalSystem compile(
        const ExecutableEquationSystem& equations,
        CompiledExecutionProgram& program,
        const FDM::NumericalRecipeSet& recipes,
        const TermProviderCatalog& catalog,
        const CompiledTimeRecipe& time,
        const std::vector<FormulaOperatorBinding>& bindings) {
    CompiledNumericalSystem result;
    result.recipes=recipes;
    if (time.id()!=recipes.time.id())
        throw std::runtime_error("Compiled TemporalMethod differs from selected time recipe.");
    result.time.recipe=time;
    for (const auto& call:program.steps) {
        result.requiredHaloWidth=std::max(result.requiredHaloWidth,call.requiredHaloWidth);
        for (const auto& workspace:call.numericalWorkspaces)
            appendUnique(result.workspaceRequirements,workspace);
    }
    const auto bindRecipe=[&](const FDM::TermRecipe& recipe,
                              RecipeConsumerKind kind,std::string consumer) {
        result.recipeBindings.push_back({recipe,kind,std::move(consumer)});
    };
    for (auto& call:program.steps)
        if (call.spatialTerms || !call.sourceMathInputs.empty()) {
            const auto first=result.operators.size();
            auto selected=bindings;
            for(const auto& op:call.numericalOperators) {
                selected.erase(std::remove_if(selected.begin(),selected.end(),[&](const auto& other){return other.formula==op.formula && other.occurrence==op.occurrence && other.operation==op.operation;}),selected.end());
                selected.push_back(op);
            }
            if(!call.lowerSpatial)
                throw std::runtime_error("Selected method has no spatial lowering implementation: "+
                    call.equationMethod+" at "+call.source.occurrence);
            try {
                call.lowerSpatial(equations,call,recipes,catalog,selected,result);
            } catch(const std::exception& error) {
                throw std::runtime_error("Spatial lowering: occurrence="+call.source.occurrence
                    +" equation="+call.source.equation+" method="+call.equationMethod+": "+error.what());
            }
            for(auto n=first;n<result.operators.size();++n)result.operators[n].callOccurrence=call.source.occurrence;
        }
    for (const auto& binding:bindings) {
        const bool addressed=std::any_of(result.operators.begin(),
            result.operators.end(),[&](const CompiledSpatialBinding& item) {
                return (binding.formula.empty()
                    && binding.operation && *binding.operation==item.operation)
                    || (binding.formula==item.formulaId
                        && (binding.occurrence.empty()
                            || binding.occurrence==item.occurrence));
            });
        if (!addressed)
            throw std::runtime_error("Numerical binding '"+binding.formula+"/"
                +binding.occurrence+"' addresses no selected Equation occurrence.");
    }
    for (const auto& operation:equations.operations)
        for (const auto role:operation.consumedRecipes) {
            const FDM::TermRecipe* recipe=role==OperationRecipeRole::Convection
                ? (recipes.convection ? &*recipes.convection : nullptr)
                : (recipes.diffusion ? &*recipes.diffusion : nullptr);
            if (!recipe)
                throw std::runtime_error("Executable operation '"
                    +operation.operation+"' consumes an unselected recipe.");
            bindRecipe(*recipe,RecipeConsumerKind::ExecutableOperation,
                       operation.operation);
        }
    if (recipes.diffusion && std::none_of(result.recipeBindings.begin(),
            result.recipeBindings.end(),[](const CompiledRecipeBinding& binding) {
                return binding.recipe.role()==FDM::TermRole::Diffusion;
            }))
        throw std::runtime_error("A diffusion TermRecipe was selected without a consuming Equation occurrence.");
    return result;
}


const FDM::TermRecipe* findUniqueRecipe(const CompiledNumericalSystem& system,
                                       FDM::TermRole role) {
    const FDM::TermRecipe* found=nullptr;
    for (const auto& binding:system.operators) {
        if (!binding.recipe || binding.recipe->role()!=role) continue;
        if (found && found->id()!=binding.recipe->id())
            throw std::runtime_error("Compiled numerical system has multiple recipes for one fused role.");
        found=&*binding.recipe;
    }
    return found;
}

const FDM::TermRecipe& requireUniqueRecipe(const CompiledNumericalSystem& system,
                                           FDM::TermRole role) {
    if (const auto* recipe=findUniqueRecipe(system,role)) return *recipe;
    throw std::runtime_error(std::string("Compiled numerical system has no bound ")
                             +FDM::toString(role)+" recipe.");
}

} // namespace SF::System::NumericalCompiler
