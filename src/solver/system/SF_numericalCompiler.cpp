/// @file SF_numericalCompiler.cpp
/// @brief Bind selected Equation AST occurrences to spatial numerical providers.

#include "SF_numericalCompiler.h"
#include "SF_config.h"

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

CompiledNumericalSystem compile(
        const ExecutableEquationSystem& equations,
        const CompiledExecutionProgram& program,
        const FDM::NumericalRecipeSet& recipes,
        const TermProviderCatalog& catalog,
        const CompiledTimeRecipe& time,
        const ProviderRegistry& providers,
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
    for (const auto& call:program.steps)
        if (call.spatialTerms || !call.sourceMathInputs.empty()) {
            auto selected=bindings;
            for(const auto& op:call.numericalOperators) {
                selected.erase(std::remove_if(selected.begin(),selected.end(),[&](const auto& other){return other.formula==op.formula && other.occurrence==op.occurrence && other.operation==op.operation;}),selected.end());
                selected.push_back(op);
            }
            if(!lowerDirectSpatial(equations,call,recipes,catalog,selected,result))
                providers.at(call.equationMethod).lowerSpatial(equations,call,recipes,catalog,selected,result);
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
