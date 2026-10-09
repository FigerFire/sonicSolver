/// @file SF_numericalCompiler.cpp
/// @brief Bind selected Equation AST occurrences to spatial numerical providers.

#include "SF_numericalCompiler.h"
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

struct SelectedFormula {
    std::string formula;
    std::string execution;
    std::string output;
    std::string method;
    bool primitiveSourceRequired = false;
    std::string residualWorkspace;
};

std::vector<SelectedFormula> selectFormulas(
        const CompiledExecutionProgram& program) {
    std::vector<SelectedFormula> selected;
    for (const auto& step:program.steps) {
        if (step.spatialTerms) {
            for (const auto& call:step.calls)
                selected.push_back({call.equation,step.source.equation,
                                    call.target,step.equationMethod,step.primitiveSourceRequired,step.residualWorkspace});
        } else if (!step.sourceMathInputs.empty()) {
            for (const auto& formula:step.sourceMathInputs)
                selected.push_back({formula,step.source.equation,
                    step.source.target.symbol,step.equationMethod,step.primitiveSourceRequired,step.residualWorkspace});
        }
        // Other methods own their local numerical realization. The spatial
        // compiler consumes the method's declared residual/source input set.
    }
    if (!program.root.children.empty() || !program.steps.empty()) return selected;
    for (const auto& call:program.legacySpatialInputs)
        selected.push_back({call.equation,call.equation,call.target,"",false,"residual"});
    return selected;
}

std::string symbolAt(const FormulaExpr& expression,std::size_t index,
                     const Equation& formula,const std::string& path) {
    if (expression.arguments.size()<=index
        || expression.arguments[index].kind!=FormulaExpr::Kind::Symbol)
        throw std::runtime_error("Equation '"+formula.id+"' occurrence '"+path
            +"' has an unsupported non-symbol operator operand.");
    return expression.arguments[index].name;
}

} // namespace

CompiledNumericalSystem compile(
        const ExecutableEquationSystem& equations,
        const CompiledExecutionProgram& program,
        const FDM::NumericalRecipeSet& recipes,
        const TermProviderCatalog& catalog,
        const CompiledTimeRecipe& time,
        const std::vector<FormulaOperatorBinding>& bindings) {
    if (recipes.time.topology()!=FDM::TimeTopology::ExplicitStages)
        throw std::runtime_error("Current Equation numerical recipes require an explicit time topology.");
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
    const auto selectedFormulas=selectFormulas(program);
    const auto containsOperator=[](const auto& self,const FormulaExpr& expression,std::string_view name)->bool {
        if (expression.kind==FormulaExpr::Kind::Operator && expression.name==name) return true;
        return std::any_of(expression.arguments.begin(),expression.arguments.end(),
            [&](const auto& child) { return self(self,child,name); });
    };
    if (recipes.diffusion) {
        const bool equationConsumer=std::any_of(selectedFormulas.begin(),selectedFormulas.end(),[&](const auto& selected) {
            const auto& equation=equations.registry.at(selected.formula);
            return containsOperator(containsOperator,equation.lhs,"diffusion")
                || containsOperator(containsOperator,equation.rhs,"diffusion");
        });
        const bool operationConsumer=std::any_of(equations.operations.begin(),equations.operations.end(),[](const auto& operation) {
            return std::find(operation.consumedRecipes.begin(),operation.consumedRecipes.end(),OperationRecipeRole::Diffusion)
                !=operation.consumedRecipes.end();
        });
        if (!equationConsumer && !operationConsumer)
            throw std::runtime_error("A diffusion TermRecipe was selected without a consuming Equation occurrence.");
    }
    const auto bindRecipe=[&](const FDM::TermRecipe& recipe,
                              RecipeConsumerKind kind,std::string consumer) {
        result.recipeBindings.push_back({recipe,kind,std::move(consumer)});
    };
    for (const auto& selected:selectedFormulas) {
        const auto& formula=equations.registry.at(selected.formula);
        int transientCount=0;
        const auto visit=[&](const auto& self,const FormulaExpr& expression,
                             CompiledSpatialBinding::Side side,
                             const std::string& path)->void {
            if (expression.kind==FormulaExpr::Kind::Subtract
                || expression.kind==FormulaExpr::Kind::Multiply
                || expression.kind==FormulaExpr::Kind::Divide
                || expression.kind==FormulaExpr::Kind::Negate)
                throw std::runtime_error("Equation '"+formula.id
                    +"' at '"+path+"' requires unsupported signed/coefficient "
                    "lowering in the fused numerical backend.");
            const bool legacySource=!formula.authored
                && side==CompiledSpatialBinding::Side::Right
                && expression.kind==FormulaExpr::Kind::Symbol;
            const bool source=expression.kind==FormulaExpr::Kind::Operator
                && expression.name=="source";
            const bool convection=expression.kind==FormulaExpr::Kind::Operator
                && expression.name=="div";
            const bool diffusion=expression.kind==FormulaExpr::Kind::Operator
                && expression.name=="diffusion";
            if ((source || legacySource)
                    && side!=CompiledSpatialBinding::Side::Right
                || (convection || diffusion)
                    && side!=CompiledSpatialBinding::Side::Left)
                throw std::runtime_error("Equation '"+formula.id+"' at '"+path
                    +"' requires unsupported term-side/sign lowering.");
            if (source || convection || diffusion || legacySource) {
                FormulaExpr sourceExpression;
                const FormulaExpr* numerical=&expression;
                if (legacySource) {
                    sourceExpression=FormulaExpr::op("source",
                        {FormulaExpr::symbol(expression.name)},path);
                    numerical=&sourceExpression;
                }
                const std::string occurrence=expression.occurrence.empty()
                    ? path : expression.occurrence;
                const std::string opName=numerical->name;
                const std::string primary=symbolAt(*numerical,
                    diffusion?1:0,formula,occurrence);
                const std::string secondary=diffusion
                    ? symbolAt(*numerical,0,formula,occurrence) : "";
                const FDM::TermRecipe* requested=nullptr;
                if (convection) {
                    if (!recipes.convection)
                        throw std::runtime_error("Equation '"+formula.id
                            +"' occurrence '"+occurrence
                            +"' has no convection TermRecipe.");
                    requested=&*recipes.convection;
                } else if (diffusion) {
                    if (!recipes.diffusion)
                        throw std::runtime_error("Equation '"+formula.id
                            +"' occurrence '"+occurrence
                            +"' has no diffusion TermRecipe.");
                    requested=&*recipes.diffusion;
                }
                const auto provider=catalog.resolve({formula,*numerical,
                    occurrence,selected.output,selected.method,requested},bindings);
                const bool isSource=source || legacySource;
                const bool primitiveSource=isSource && selected.primitiveSourceRequired
                    && provider.primitiveSource;
                if (provider.status!=BindingStatus::Resolved
                    || (!provider.recipe && !primitiveSource && !provider.stageSource.evaluate))
                    throw std::runtime_error("Equation provider "
                        +std::string(provider.status==BindingStatus::Invalid
                            ? "Invalid: " : "Unsupported: ")
                        +(provider.reason.empty()
                            ? "missing recipe for '"+occurrence+"'"
                            : provider.reason));
                if (isSource && !selected.primitiveSourceRequired
                    && !provider.conservativeSource && !provider.stageSource.evaluate)
                    throw std::runtime_error("Conservative source provider '"
                        +provider.id+"' has no compiled execution kernel.");
                result.operators.emplace_back(formula.id,selected.execution,
                    selected.output,selected.method,occurrence,opName,
                    primary,secondary,side,provider.recipe,provider.id,
                    provider.owner,provider.compiledDataAvailable,
                    provider.primitiveSource,provider.conservativeSource,
                    provider.primitiveSpatial);
                result.operators.back().stageSource=provider.stageSource;
                appendUnique(result.providerRequirements,provider.id);
                if (requested) {
                    bindRecipe(*requested,RecipeConsumerKind::EquationTerm,
                        formula.id+"@"+occurrence);
                    if (convection) {
                        result.requiredHaloWidth=std::max(
                            result.requiredHaloWidth,requested->haloWidth());
                        appendUnique(result.workspaceRequirements,selected.residualWorkspace);
                        if (requested->reconstruction()
                            ==FDM::ReconstructionVariable::Characteristic)
                            appendUnique(result.workspaceRequirements,
                                         "characteristicReconstruction");
                    } else {
                        appendUnique(result.workspaceRequirements,
                            selected.residualWorkspace);
                    }
                } else {
                    if (provider.recipe)
                        bindRecipe(*provider.recipe,RecipeConsumerKind::EquationTerm,
                            formula.id+"@"+occurrence);
                    appendUnique(result.workspaceRequirements,
                        selected.residualWorkspace);
                }
                return;
            }
            if (expression.kind==FormulaExpr::Kind::Operator
                && expression.name!="ddt")
                throw std::runtime_error("Equation '"+formula.id
                    +"' has unsupported numerical operator '"+expression.name
                    +"' at '"+path+"'.");
            if (expression.kind==FormulaExpr::Kind::Operator
                && expression.name=="ddt") {
                if (++transientCount!=1 || expression.arguments.size()!=1
                    || side!=CompiledSpatialBinding::Side::Left
                    || symbolAt(expression,0,formula,path)!=selected.output)
                    throw std::runtime_error("Equation ddt does not match HOW Output.");
                return;
            }
            if (expression.kind==FormulaExpr::Kind::Symbol
                || (expression.kind==FormulaExpr::Kind::Constant
                    && expression.constant!=0.0))
                throw std::runtime_error("Equation '"+formula.id+"' at '"+path
                    +"' has unsupported unbound residual data.");
            for (std::size_t i=0;i<expression.arguments.size();++i)
                self(self,expression.arguments[i],side,
                    path+"/"+std::to_string(i));
        };
        visit(visit,formula.lhs,CompiledSpatialBinding::Side::Left,"lhs");
        visit(visit,formula.rhs,CompiledSpatialBinding::Side::Right,"rhs");
        if (transientCount!=1)
            throw std::runtime_error("Fused transport Equation requires exactly one ddt(Output).");
    }
    for (const auto& binding:bindings) {
        const bool addressed=std::any_of(result.operators.begin(),
            result.operators.end(),[&](const CompiledSpatialBinding& item) {
                return (binding.formula.empty()
                    && binding.occurrence==item.operatorName)
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
        const NumericalSelection& selected,const std::vector<LegacyPlanFragment>& legacyFragments,
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
    compiledExecution.legacySpatialInputs=selected.legacySpatialInputs;
    TermProviderCatalog termProviders=TermProviderCatalog::builtIn();
    for (const auto& descriptor:selected.termProviders) termProviders.addSource(descriptor);
    result.numericalSystem=compile(result.executableSystem,compiledExecution,
        selected.recipes,termProviders,compiledTime,selected.operatorBindings);
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
        result.executableSystem,result.executionPolicies,
        result.numericalSystem.time.recipe,
        legacyFragments,compiledExecution);
    result.solvePlan.sourceProgram=executionProgram;
    result.solvePlan.sourceProgram.root=compiledExecution.root;
    result.solvePlan.compiledProgram=std::move(compiledExecution);
    result.runtime.capabilities=compileExecutionCapabilities(
        result.executableSystem,result.executionPolicies,&result.solvePlan);
    result.runtime.operationBindings=compileOperationBindings(result.executableSystem,
        result.numericalSystem,result.solvePlan,result.executionPolicies,additionalContributions);
    result.runtime.report=reportOperationBindings(result.solvePlan,result.runtime.operationBindings);
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
