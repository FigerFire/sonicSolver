/// @file SF_numericalCompiler.cpp
/// @brief Bind selected Equation AST occurrences to spatial numerical providers.

#include "SF_numericalCompiler.h"
#include "SF_formulaCompiler.h"
#include "SF_executionContract.h"
#include "SF_methodObjects.h"
#include "SF_solvePlan.h"
#include "SF_providerResolver.h"
#include "SF_stateRealization.h"

#include <algorithm>
#include <functional>
#include <stdexcept>
#include <tuple>

namespace SF::System::NumericalCompiler {
namespace {

void appendUnique(std::vector<std::string>& values,std::string value) {
    if (std::find(values.begin(),values.end(),value)==values.end())
        values.push_back(std::move(value));
}

} // namespace

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
        if (call.spatialTerms || !call.sourceMathInputs.empty())
            providers.at(call.equationMethod).lowerSpatial(equations,call,recipes,catalog,bindings,result);
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

void compileSystem(ResolvedSimulationSystem& result,const ExecutionProgram& executionProgram,
        const NumericalSelection& selected,
        bool additionalContributions) {
    result.realization=compileStateRealization(result.executableSystem.state,result.executableSystem.constraints);
    const auto temporalMethods=builtinTemporalMethods();
    const auto& temporalMethod=temporalMethods.at(selected.recipes.time.id());
    const auto compiledTime=temporalMethod.compile(selected.recipes.time);
    auto compiledExecution=compileExecutionProgram(result.executableSystem,result.executableSystem.state,executionProgram,
        selected.bindings,builtinProviders(),&compiledTime,&temporalMethod);
    for (const auto& id:result.executableSystem.state.solutionVariables()) {
        const auto& symbol=result.executableSystem.state.at(id);
        const bool owned=std::any_of(compiledExecution.steps.begin(),compiledExecution.steps.end(),[&](const auto& step) {
            if (std::find(step.writes.begin(),step.writes.end(),id)!=step.writes.end()) return true;
            return std::any_of(step.target.resources.begin(),step.target.resources.end(),[&](const auto& resource) {
                return resource.storage==symbol.storageKey && resource.access!=ResourceAccessMode::Read
                    && resource.componentOffset<=symbol.componentOffset
                    && resource.componentOffset+resource.components>=symbol.componentOffset+symbol.components;
            });
        });
        if (!owned) throw std::runtime_error("No compiled HOW output owns selected solution STATE '"+id+"'.");
    }
    TermProviderCatalog termProviders=TermProviderCatalog::builtIn();
    for (const auto& descriptor:selected.termProviders) termProviders.addSource(descriptor);
    result.numericalSystem=compile(result.executableSystem,compiledExecution,
        selected.recipes,termProviders,compiledTime,builtinProviders(),selected.operatorBindings);
    result.numericalSystem.dt=selected.dt;
    result.numericalSystem.phaseTransport=selected.phaseTransport;
    result.numericalSystem.pressureFaceCoupling=selected.pressureFaceCoupling;
    result.numericalSystem.pressureOperator=selected.pressureOperator;
    // Compiled operation/resource declarations originate in selected providers.
    for (const auto& call:compiledExecution.steps) {
        for (const auto& operation:call.operations) {
            const auto duplicate=std::find_if(result.executableSystem.operations.begin(),result.executableSystem.operations.end(),
                [&](const auto& existing) { return existing.operation==operation.operation; });
            if (duplicate==result.executableSystem.operations.end()) result.executableSystem.operations.push_back(operation);
        }
    }
    for (auto& step:compiledExecution.steps) {
        for (const auto& call:step.calls) {
            for (const auto& term:result.numericalSystem.operators) {
                if ((term.formulaId!=call.equation
                        && std::find(step.sourceMathInputs.begin(),
                            step.sourceMathInputs.end(),term.formulaId)
                           ==step.sourceMathInputs.end())
                    || term.provider.empty()) continue;
                if (std::find(step.operatorBindings.begin(),
                        step.operatorBindings.end(),term.provider)
                    ==step.operatorBindings.end())
                    step.operatorBindings.push_back(term.provider);
            }
        }
        if (step.temporalResidual && step.spatialTerms) {
            if (step.spatialTerms && step.operatorBindings.empty())
                throw std::runtime_error("ConservativeResidual has no bound spatial provider.");
            step.requirements.push_back("halo depth="
                +std::to_string(result.numericalSystem.requiredHaloWidth));
            step.requirements.push_back("canonical face COPY before residual assembly");
        }
    }
    result.solvePlan = SolvePlanner::compile(
        result.executableSystem,
        result.numericalSystem.time.recipe,
        compiledExecution);
    result.solvePlan.sourceProgram=executionProgram;
    result.solvePlan.sourceProgram.root=compiledExecution.root;
    result.solvePlan.compiledProgram=std::move(compiledExecution);
    result.runtime.capabilities=compileExecutionCapabilities(
        result.executableSystem,&result.solvePlan);
    result.runtime.operationBindings=compileOperationBindings(result.executableSystem,
        result.numericalSystem,result.solvePlan,additionalContributions);
    result.runtime.report=reportOperationBindings(result.solvePlan,result.runtime.operationBindings);
    validateMathematicalCoverage(result.executableSystem,result.solvePlan,result.runtime.report.status==RuntimeStatus::Runnable);
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

std::vector<FDM::SourceKind> sourceKinds(const CompiledNumericalSystem& system) {
    std::vector<FDM::SourceKind> result;
    for (const auto& binding:system.operators) {
        if (!binding.recipe
            || binding.recipe->role()!=FDM::TermRole::Source) continue;
        const auto source=binding.recipe->source();
        if (std::find(result.begin(),result.end(),source)==result.end())
            result.push_back(source);
    }
    return result;
}

} // namespace SF::System::NumericalCompiler
