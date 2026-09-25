/// @file SF_numericalCompiler.cpp
/// @brief Mathematical term 到 numerical recipe 的 production lowering。

#include "SF_numericalCompiler.h"

#include "SF_config.h"

#include <algorithm>
#include <stdexcept>

namespace SF::System::NumericalCompiler {
namespace {

bool isCoreFluidEquation(const std::string& id) {
    return id == "E_MASS" || id == "E_MOMENTUM" || id == "E_ENERGY";
}

void appendUnique(std::vector<std::string>& values, std::string value) {
    if (std::find(values.begin(),values.end(),value) == values.end()) {
        values.push_back(std::move(value));
    }
}

FDM::TermRecipe sourceRecipe(const std::string& symbol) {
    if (symbol == "gravity") {
        return FDM::builtInSourceRecipe(FDM::SourceKind::Gravity);
    }
    if (symbol == "MRF") {
        return FDM::builtInSourceRecipe(FDM::SourceKind::MRF);
    }
    if (symbol == "wallHeat") {
        return FDM::builtInSourceRecipe(FDM::SourceKind::WallHeat);
    }
    throw std::runtime_error(
        "No built-in source TermRecipe is registered for mathematical source '"
        +symbol+"'.");
}

} // namespace

CompiledNumericalSystem compile(
        const ExecutableEquationSystem& equations,
        const FDM::NumericalRecipeSet& recipes) {
    if (recipes.time.topology() != FDM::TimeTopology::ExplicitStages) {
        throw std::runtime_error(
            "Current TermRecipes provide ExplicitResidual only; selected "
            "TimeRecipe has an unsupported topology.");
    }

    CompiledNumericalSystem result;
    result.recipes = recipes;
    // 编译后的时间策略与 term recipe 来自同一份输入，因此 runtime 不需要
    // 再读 raw numerics 选择。
    result.time.recipe = recipes.time;
    const auto bindRecipe = [&](const FDM::TermRecipe& recipe,
                                RecipeConsumerKind kind,
                                std::string consumer) {
        result.recipeBindings.push_back({recipe,kind,std::move(consumer)});
    };

    for (const auto& equation : equations.equations) {
        if (!isCoreFluidEquation(equation.id)) continue;
        const auto& definition = equations.equationDefinitions.at(equation.id);
        std::size_t ordinal = 0;
        const auto bind = [&](const Equation::Expression& expression) {
            for (const auto& term : expression.terms) {
                const std::size_t currentOrdinal = ordinal++;
                if (term.kind == Equation::TermKind::Divergence) {
                    if (!recipes.convection) {
                        throw std::runtime_error(
                            "Equation '"+equation.id+"' contains convection term '"
                            +term.primary.name
                            +"' but no convection TermRecipe was selected.");
                    }
                    result.terms.emplace_back(
                        equation.id,term,currentOrdinal,*recipes.convection);
                    bindRecipe(*recipes.convection,RecipeConsumerKind::EquationTerm,
                               equation.id+"["+std::to_string(currentOrdinal)+"]");
                    result.requiredHaloWidth = std::max(
                        result.requiredHaloWidth,
                        recipes.convection->haloWidth());
                    appendUnique(result.workspaceRequirements,"faceFlux");
                    appendUnique(result.workspaceRequirements,"residual");
                    appendUnique(
                        result.workspaceRequirements,"characteristicReconstruction");
                    appendUnique(
                        result.providerRequirements,
                        std::string("term.convection.")
                            +FDM::toString(recipes.convection->id()));
                } else if (term.kind == Equation::TermKind::Diffusion) {
                    if (!recipes.diffusion) {
                        throw std::runtime_error(
                            "Equation '"+equation.id+"' contains diffusion term '"
                            +term.primary.name
                            +"' but no diffusion TermRecipe was selected.");
                    }
                    result.terms.emplace_back(
                        equation.id,term,currentOrdinal,*recipes.diffusion);
                    bindRecipe(*recipes.diffusion,RecipeConsumerKind::EquationTerm,
                               equation.id+"["+std::to_string(currentOrdinal)+"]");
                    appendUnique(result.workspaceRequirements,"residual");
                    appendUnique(
                        result.providerRequirements,
                        std::string("term.diffusion.")
                            +FDM::toString(recipes.diffusion->id()));
                } else if (term.kind == Equation::TermKind::Source
                           && term.primary.name != "zero") {
                    const auto recipe = sourceRecipe(term.primary.name);
                    result.terms.emplace_back(
                        equation.id,term,currentOrdinal,recipe);
                    bindRecipe(recipe,RecipeConsumerKind::EquationTerm,
                               equation.id+"["+std::to_string(currentOrdinal)+"]");
                    appendUnique(result.workspaceRequirements,"residual");
                    appendUnique(
                        result.providerRequirements,
                        std::string("term.source.")+FDM::toString(recipe.id()));
                }
            }
        };
        bind(definition.left);
        bind(definition.right);
    }

    for (const ExecutableOperation& operation : equations.operations) {
        for (OperationRecipeRole role : operation.consumedRecipes) {
            const FDM::TermRecipe* recipe = role == OperationRecipeRole::Convection
                ? (recipes.convection ? &*recipes.convection : nullptr)
                : (recipes.diffusion ? &*recipes.diffusion : nullptr);
            if (!recipe) {
                throw std::runtime_error(
                    "Executable operation '"+operation.operation
                    +"' consumes a numerical recipe that was not selected.");
            }
            bindRecipe(*recipe,RecipeConsumerKind::ExecutableOperation,
                       operation.operation);
        }
    }
    const bool usedDiffusion = std::any_of(
        result.recipeBindings.begin(),result.recipeBindings.end(),
        [](const CompiledRecipeBinding& binding) {
            return binding.recipe.role() == FDM::TermRole::Diffusion;
        });
    if (recipes.diffusion && !usedDiffusion) {
        throw std::runtime_error(
            "A diffusion TermRecipe was selected but the Equation System "
            "contains no diffusion term or consuming operation.");
    }
    return result;
}

const FDM::TermRecipe* findUniqueRecipe(
        const CompiledNumericalSystem& system, FDM::TermRole role) {
    const FDM::TermRecipe* found = nullptr;
    for (const auto& term : system.terms) {
        if (term.recipe.role() != role) continue;
        if (found && found->id() != term.recipe.id()) {
            throw std::runtime_error(
                "Compiled numerical system contains multiple recipes for one "
                "fused term role.");
        }
        found = &term.recipe;
    }
    return found;
}

const FDM::TermRecipe& requireUniqueRecipe(
        const CompiledNumericalSystem& system, FDM::TermRole role) {
    const auto* found = findUniqueRecipe(system,role);
    if (found) return *found;
    throw std::runtime_error(
        std::string("Compiled numerical system has no bound ")
        +FDM::toString(role)+" recipe.");
}

std::vector<FDM::SourceKind> sourceKinds(
        const CompiledNumericalSystem& system) {
    std::vector<FDM::SourceKind> result;
    for (const auto& term : system.terms) {
        if (term.recipe.role() != FDM::TermRole::Source) continue;
        const auto source = term.recipe.source();
        if (std::find(result.begin(),result.end(),source) == result.end()) {
            result.push_back(source);
        }
    }
    return result;
}

} // namespace SF::System::NumericalCompiler
