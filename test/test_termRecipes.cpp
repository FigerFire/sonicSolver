#include "SF_config.h"
#include "solver/system/SF_numericalCompiler.h"
#include "solver/system/SF_systemBuilder.h"
#include "solver/system/SF_systemPrinter.h"
#include "models/physics/SF_sourceContribution.h"
#include "models/physics/mrf/SF_frameProvider.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include "app/application/model/SF_configParser.h"

namespace {

using namespace SF;

[[noreturn]] void fail(const std::string& message) {
    std::cerr << "Term recipe test failed: " << message << '\n';
    std::exit(1);
}

void require(bool condition, const std::string& message) {
    if (!condition) fail(message);
}

std::string rawSignature(const System::RawEquationSystem& raw) {
    std::ostringstream output;
    for (const auto& equation:raw.registry.entries()) output << System::formulaText(equation) << '\n';
    return output.str();
}

bool hasTerm(const System::RawEquationSystem& raw,
             const std::string& equationId,
             SF::Equation::TermKind kind) {
    const auto& equation=raw.registry.at(equationId);
    const auto visit=[&](const auto& self,const System::FormulaExpr& expression)->bool {
        if (kind==SF::Equation::TermKind::Diffusion && expression.kind==System::FormulaExpr::Kind::Operator
            && expression.name=="diffusion") return true;
        return std::any_of(expression.arguments.begin(),expression.arguments.end(),[&](const auto& child) { return self(self,child); });
    };
    return visit(visit,equation.lhs) || visit(visit,equation.rhs);
}

FDM::SolverConfig explicitConfig(FDM::TermRecipe convection) {
    FDM::SolverConfig config;
    config.numerics.timeRecipe =
        FDM::builtInTimeRecipe(FDM::TimeRecipeId::ForwardEuler);
    config.numerics.recipes.time = config.numerics.timeRecipe;
    config.numerics.recipes.convection = convection;
    config.numerics.recipes.diffusion.reset();
    config.numerics.termRecipesDeclared = true;
    return config;
}

template<class Callable>
bool throwsContaining(Callable&& callable, const std::string& expected) {
    try {
        callable();
    } catch (const std::exception& error) {
        const bool matched=std::string(error.what()).find(expected)!=std::string::npos;
        if (!matched) std::cerr << "Expected error: " << expected << "; actual: " << error.what() << '\n';
        return matched;
    }
    return false;
}

} // namespace

int main() {
    System::BuildRequest inviscid;
    inviscid.composition.stateDeclared = true;
    inviscid.composition.solutionVariables = {"rho","rhoU","rhoE"};
    inviscid.templateOrigin = System::PhysicsTemplateKind::SingleFluid;
    inviscid.singleFluidPreset = System::SingleFluidPresetSpec{false};

    const auto rusanovRecipe =
        FDM::resolveConvectionTermRecipe("weno7Rusanov");
    const auto stegerRecipe =
        FDM::resolveConvectionTermRecipe("weno7Steger");
    const auto rusanov = System::build(
        explicitConfig(rusanovRecipe),inviscid);
    auto overrideRequest=inviscid;
    overrideRequest.operatorBindings={{"","div","convection.conservativeFlux"},
        {"momentum","","convection.conservativeFlux"},
        {"momentum","momentum.convection","convection.conservativeFlux"}};
    const auto overridden=System::build(explicitConfig(rusanovRecipe),overrideRequest);
    require(overridden.numericalSystem.operators.size()
                ==rusanov.numericalSystem.operators.size(),
            "Equation override changed the physical operator collection.");
    overrideRequest.operatorBindings.back().provider="missing.provider";
    require(throwsContaining([&] {
        (void)System::build(explicitConfig(rusanovRecipe),overrideRequest);
    },"missing.provider"),"Explicit provider selection silently fell back.");
    overrideRequest.operatorBindings={{"momentum","absent.occurrence",
                                       "convection.conservativeFlux"}};
    require(throwsContaining([&] {
        (void)System::build(explicitConfig(rusanovRecipe),overrideRequest);
    },"addresses no selected Equation occurrence"),
        "Unused Equation occurrence binding was silently ignored.");
    auto signedSystem=rusanov.executableSystem;
    auto signedFormula=signedSystem.registry.at("momentum");
    signedFormula.lhs=System::FormulaExpr::negate(signedFormula.lhs);
    signedSystem.registry.replace(signedFormula);
    require(throwsContaining([&] {
        (void)System::NumericalCompiler::compile(signedSystem,
            rusanov.solvePlan.compiledProgram,
            explicitConfig(rusanovRecipe).numerics.recipes,
            System::TermProviderCatalog::builtIn(),rusanov.numericalSystem.time.recipe);
    },"unsupported signed/coefficient"),
        "Fused backend silently ignored a Equation sign.");
    auto unknownFluxSystem=rusanov.executableSystem;
    auto unknownFlux=unknownFluxSystem.registry.at("continuity");
    unknownFlux.lhs.arguments.back().arguments.front()=
        System::FormulaExpr::symbol("userFluxWithoutProvider");
    unknownFluxSystem.registry.replace(unknownFlux);
    require(throwsContaining([&] {
        (void)System::NumericalCompiler::compile(unknownFluxSystem,
            rusanov.solvePlan.compiledProgram,
            explicitConfig(rusanovRecipe).numerics.recipes,
            System::TermProviderCatalog::builtIn(),rusanov.numericalSystem.time.recipe);
    },"Unsupported"),"Fused provider executed an unknown mathematical flux.");
    auto duplicateTimeSystem=rusanov.executableSystem;
    auto duplicateTime=duplicateTimeSystem.registry.at("continuity");
    duplicateTime.lhs=System::FormulaExpr::add(duplicateTime.lhs,
        System::FormulaExpr::op("ddt",{System::FormulaExpr::symbol("rho")}));
    duplicateTimeSystem.registry.replace(duplicateTime);
    require(throwsContaining([&] {
        (void)System::NumericalCompiler::compile(duplicateTimeSystem,
            rusanov.solvePlan.compiledProgram,
            explicitConfig(rusanovRecipe).numerics.recipes,
            System::TermProviderCatalog::builtIn(),rusanov.numericalSystem.time.recipe);
    },"ddt"),"Fused transport silently ignored a duplicate time derivative.");
    const auto steger = System::build(
        explicitConfig(stegerRecipe),inviscid);

    require(rawSignature(rusanov.rawSystem) == rawSignature(steger.rawSystem),
            "spatial recipe changed RawEquationSystem physics");
    require(std::any_of(
                rusanov.rawSystem.contributions.begin(),
                rusanov.rawSystem.contributions.end(),
                [](const System::ContributionRecord& contribution) {
                    return contribution.id == "builtin.singleFluidNavierStokes";
                }),
            "density core equations were not installed by the single-fluid preset");
    require(!hasTerm(rusanov.rawSystem,"momentum",
                     SF::Equation::TermKind::Diffusion)
                && !hasTerm(rusanov.rawSystem,"energy",
                            SF::Equation::TermKind::Diffusion),
            "inviscid preset contains a physical diffusion term");
    require(rusanov.timeRecipe.id() == steger.timeRecipe.id()
                && rusanov.timeRecipe.stageCount()
                    == steger.timeRecipe.stageCount(),
            "spatial recipe changed TimeRecipe topology");
    require(rusanov.numericalSystem.requiredHaloWidth == 4,
            "WENO7 recipe did not own its four-cell halo requirement");
    require(System::NumericalCompiler::requireUniqueRecipe(
                rusanov.numericalSystem,FDM::TermRole::Convection).id()
                == FDM::TermRecipeId::Weno7Rusanov,
            "Rusanov recipe was not bound into CompiledNumericalSystem");
    require(System::NumericalCompiler::requireUniqueRecipe(
                steger.numericalSystem,FDM::TermRole::Convection).id()
                == FDM::TermRecipeId::Weno7Steger,
            "Steger recipe was not bound into CompiledNumericalSystem");
    require(rusanov.numericalSystem.operators.size() == 3,
            "core mass/momentum/energy convection terms were not all bound");
    require(System::describe(rusanov).find(
                "momentum@momentum.convection -> rhoU div(momentumFlux) -> weno7Rusanov")
                != std::string::npos,
            "explain output does not expose compiled term binding");

    auto missingConvection = explicitConfig(rusanovRecipe);
    missingConvection.numerics.recipes.convection.reset();
    require(throwsContaining(
                [&] { (void)System::build(missingConvection,inviscid); },
                "no convection TermRecipe"),
            "missing convection recipe did not fail fast");

    auto unusedDiffusion = explicitConfig(rusanovRecipe);
    unusedDiffusion.numerics.recipes.diffusion =
        FDM::builtInDiffusionRecipe(FDM::ViscousScheme::Central2);
    require(throwsContaining(
                [&] { (void)System::build(unusedDiffusion,inviscid); },
                "without a consuming Equation occurrence"),
            "unused diffusion recipe was silently accepted");

    auto viscousRequest = inviscid;
    viscousRequest.singleFluidPreset = System::SingleFluidPresetSpec{true};
    require(throwsContaining(
                [&] { (void)System::build(
                    explicitConfig(rusanovRecipe),viscousRequest); },
                "no diffusion TermRecipe"),
            "physical diffusion without a recipe did not fail fast");

    auto viscousConfig = explicitConfig(rusanovRecipe);
    viscousConfig.numerics.recipes.diffusion =
        FDM::builtInDiffusionRecipe(FDM::ViscousScheme::Central2);
    const auto viscous = System::build(viscousConfig,viscousRequest);
    require(rawSignature(viscous.rawSystem) != rawSignature(rusanov.rawSystem),
            "viscous physics did not add diffusion mathematical terms");
    require(System::NumericalCompiler::requireUniqueRecipe(
                viscous.numericalSystem,FDM::TermRole::Diffusion).id()
                == FDM::TermRecipeId::Central2Explicit,
            "physical diffusion was not bound to central2Explicit");
    require(std::any_of(viscous.numericalSystem.recipeBindings.begin(),
                        viscous.numericalSystem.recipeBindings.end(),
                        [](const System::CompiledRecipeBinding& binding) {
                            return binding.recipe.role() == FDM::TermRole::Diffusion
                                && binding.kind == System::RecipeConsumerKind::EquationTerm;
                        }),
            "diffusion term has no compiled recipe consumer");

    auto operationSystem = rusanov.executableSystem;
    System::ExecutableOperation consumingOperation;
    consumingOperation.operation = "test.diffusion.consume";
    consumingOperation.consumedRecipes = {System::OperationRecipeRole::Diffusion};
    operationSystem.operations.push_back(std::move(consumingOperation));
    auto operationRecipes = explicitConfig(rusanovRecipe).numerics.recipes;
    operationRecipes.diffusion =
        FDM::builtInDiffusionRecipe(FDM::ViscousScheme::Central2);
    const auto operationNumerics = System::NumericalCompiler::compile(
        operationSystem,rusanov.solvePlan.compiledProgram,operationRecipes,
        System::TermProviderCatalog::builtIn(),
        rusanov.numericalSystem.time.recipe);
    require(std::any_of(operationNumerics.recipeBindings.begin(),
                        operationNumerics.recipeBindings.end(),
                        [](const System::CompiledRecipeBinding& binding) {
                            return binding.kind
                                    == System::RecipeConsumerKind::ExecutableOperation
                                && binding.consumerId == "test.diffusion.consume";
                        }),
            "operation recipe consumer was not compiled");

    auto constraintRequest = inviscid;
    constraintRequest.singleFluidPreset.reset();
    constraintRequest.pressureConstraint = System::PressureConstraintSpec{false};
    require(throwsContaining(
                [&] { (void)System::build(unusedDiffusion,constraintRequest); },
                "without a consuming Equation occurrence"),
            "a constraint incorrectly exempted an unused diffusion recipe");
    require(hasTerm(viscous.rawSystem,"momentum",
                    SF::Equation::TermKind::Diffusion)
                && hasTerm(viscous.rawSystem,"energy",
                           SF::Equation::TermKind::Diffusion),
            "viscous preset did not author both physical diffusion terms");

    auto sourceConfig = explicitConfig(rusanovRecipe);
    sourceConfig.sources.enabled = {FDM::SourceKind::MRF};
    sourceConfig.sources.rotating.push_back({});
    auto sourceRequest=inviscid;
    System::SystemContribution sourceContribution;
    Physics::SourceContribution::contribute(
        sourceContribution,sourceConfig.sources);
    sourceRequest.modelContributions.push_back(std::move(sourceContribution));
    for (auto descriptor:Physics::MRF::termProviders(
             sourceConfig.sources.rotating))
        sourceRequest.termProviders.push_back(std::move(descriptor));
    const auto sourced = System::build(sourceConfig,sourceRequest);
    const auto sources =
        System::NumericalCompiler::sourceKinds(sourced.numericalSystem);
    require(sources == std::vector<FDM::SourceKind>{FDM::SourceKind::MRF},
            "MRF equation source was not compiled into one source binding");

    require(FDM::resolveDiffusionTermRecipe("central2Explicit").id()
                == FDM::TermRecipeId::Central2Explicit,
            "native diffusion recipe resolver changed identity");
    require(throwsContaining(
                [] { (void)FDM::resolveConvectionTermRecipe("notARecipe"); },
                "Unsupported convection term recipe"),
            "unknown native recipe silently fell back");

    std::cout << "Term recipe authority contract passed\n";
    return 0;
}
