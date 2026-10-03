#pragma once
#include "SF_numericalSystem.h"
#include "SF_termProviderCatalog.h"
#include "core/system/SF_numericalBinding.h"
#include "core/system/SF_solveProgram.h"
namespace SF::System { struct BuildRequest; }
namespace SF::System {
/// Composed WHICH; immutable input to binding/lowering. No execution topology.
struct NumericalSelection {
    std::vector<NumericalBinding> bindings;
    std::vector<CompiledMathRef> legacySpatialInputs;
    FDM::NumericalRecipeSet recipes;
    TimeStepPolicy dt;
    CompiledPhaseTransport phaseTransport;
    PressureFaceCoupling pressureFaceCoupling = PressureFaceCoupling::None;
    std::optional<PressureNumericalConfig> pressureOperator;
    std::vector<SourceTermProviderDescriptor> termProviders;
    std::vector<FormulaOperatorBinding> operatorBindings;
};
/// Builtin numerical defaults are composed before compilation begins.
NumericalSelection selectNumerics(const FDM::SolverConfig&, const BuildRequest&,
    const ExecutableEquationSystem&, std::vector<NumericalBinding>);
}

