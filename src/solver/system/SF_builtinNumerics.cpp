#include "SF_numericalCompiler.h"
#include "SF_buildRequest.h"
#include "SF_legacyNumerics.h"
#include "SF_methodObjects.h"
#include <algorithm>
namespace SF::System {
NumericalSelection selectNumerics(const FDM::SolverConfig& config,const BuildRequest& request,
        const ExecutableEquationSystem& equations,std::vector<NumericalBinding> bindings) {
    NumericalSelection selection;
    selection.bindings=std::move(bindings);
    if (selection.bindings.empty()) selection.legacySpatialInputs=Legacy::spatialInputs(equations);
    selection.recipes=config.numerics.recipes;
    selection.recipes.time=config.numerics.timeRecipe;
    const auto diffusion=[](const auto& self,const FormulaExpr& expression)->bool {
        if (expression.kind==FormulaExpr::Kind::Operator && expression.name=="diffusion") return true;
        return std::any_of(expression.arguments.begin(),expression.arguments.end(),
            [&](const auto& child) { return self(self,child); });
    };
    const auto providers=builtinProviders();
    const bool consumer=std::any_of(selection.bindings.begin(),selection.bindings.end(),[&](const auto& binding) {
        if (!providers.at(binding.method).usesSpatialRecipes()) return false;
        if (!equations.registry.contains(binding.equation)) return false;
        const auto& math=equations.registry.at(binding.equation);
        return diffusion(diffusion,math.lhs) || diffusion(diffusion,math.rhs);
    });
    if (consumer && !selection.recipes.diffusion && !config.numerics.termRecipesDeclared)
        selection.recipes.diffusion=FDM::builtInDiffusionRecipe(config.numerics.viscous);
    selection.dt.cfl=config.numerics.cfl;
    selection.dt.maxDeltaT=config.numerics.maxDeltaT;
    selection.phaseTransport.convection=config.pressure.phaseTransport.convection;
    selection.phaseTransport.sourceCfl=config.pressure.phaseTransport.sourceCfl;
    selection.termProviders=request.termProviders;
    selection.operatorBindings=request.operatorBindings;
    if (std::any_of(selection.bindings.begin(),selection.bindings.end(),[](const auto& binding) {
            return binding.method=="PressureMomentum";
        })) {
        selection.pressureFaceCoupling=PressureFaceCoupling::RhieChow;
        PressureNumericalConfig pressure;
        pressure.dynamicViscosity=config.numerics.dynamicViscosity;
        pressure.velocityRelaxation=config.pressure.coupling.momentumRelaxation;
        pressure.pressureRelaxation=config.pressure.coupling.pressureRelaxation;
        pressure.relativeTolerance=config.pressure.relativeTolerance;
        pressure.absoluteTolerance=config.pressure.absoluteTolerance;
        pressure.reference=config.pressure.reference;
        pressure.pressureLinear=config.pressure.linear.pressure;
        pressure.velocityBoundary=config.boundaries.velocity;
        pressure.pressureBoundary=config.boundaries.pressure;
        selection.pressureOperator=std::move(pressure);
    }
    return selection;
}
} // namespace SF::System
