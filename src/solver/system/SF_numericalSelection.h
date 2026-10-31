#pragma once
#include "SF_numericalSystem.h"
#include "SF_termProviderCatalog.h"
#include "core/system/SF_numericalBinding.h"
#include "core/system/SF_solveProgram.h"
namespace SF::System {
struct BuildRequest;
/// Composed WHICH; immutable input to binding/lowering. No execution topology.
struct NumericalSelection {
    std::vector<NumericalBinding> bindings;
    FDM::NumericalRecipeSet recipes;
    /// @brief 未显式声明离散时的既有配置默认；仅实际 spatial consumer 才会选用。
    std::optional<FDM::TermRecipe> defaultDiffusion;
    TimeStepPolicy dt;
    CompiledPhaseTransport phaseTransport;
    PressureFaceCoupling pressureFaceCoupling = PressureFaceCoupling::None;
    std::optional<PressureNumericalConfig> pressureOperator;
    std::optional<PhaseNumericalConfig> phaseOperator;
    std::vector<SourceTermProviderDescriptor> termProviders;
    std::vector<FormulaOperatorBinding> operatorBindings;
};
/// @brief 一个 occurrence 的有效 WHICH；来源独立于实际数值匹配。
struct ResolvedBinding { NumericalBinding value;std::vector<Provenance> origins; };
ResolvedBinding resolveNumerics(const EquationCall&,const std::vector<NumericalBinding>&,bool block=false,const ExecutableEquationSystem* system=nullptr);
void resolveOperatorReferences(const ExecutableEquationSystem&,const EquationCall&,NumericalBinding&);
/// Builtin numerical defaults are composed before compilation begins.
NumericalSelection selectNumerics(const FDM::SolverConfig&, const BuildRequest&,
    const ExecutableEquationSystem&, std::vector<NumericalBinding>);
}

