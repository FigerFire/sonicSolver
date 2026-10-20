/// @file SF_presets.cpp
/// @brief COMPOSE — built-in 物理方程 bundle 的四模块贡献。
///
/// 默认 equation calls 与 numerical bindings 随数学贡献声明；全局 coupling/time
/// 仍来自独立显式请求，Equation pack 不选择 loop topology。

#include "SF_systemBuilder.h"
#include "SF_equationContribution.h"
#include "SF_singleFluidPreset.h"
#include "SF_eulerianRelations.h"

#include "SF_config.h"

#include <algorithm>
#include <string>
#include <utility>
#include <tuple>
#include <vector>

namespace SF::System {
namespace Compose {
namespace {

void addState(
        SystemCompositionBuilder& system,
        std::string id,
        std::string name,
        int components = 1,
        StateRole role = StateRole::Primary,
        StorageBinding binding = StorageBinding::SpecializedExecutor,
        std::string storageKey = {},
        int componentOffset = 0,
        std::string nameSpace = {}) {
    StateSymbol unknown;
    unknown.id = std::move(id);
    unknown.name = std::move(name);
    unknown.components = components;
    unknown.shape = components == 1 ? ValueShape::Scalar : ValueShape::Vector;
    unknown.role = role;
    unknown.storageBinding = binding;
    unknown.storageKey = std::move(storageKey);
    unknown.componentOffset = componentOffset;
    unknown.nameSpace = std::move(nameSpace);
    system.addState(std::move(unknown));
}

} // namespace

void addPressureConstraintFluid(
        SystemCompositionBuilder& system,
        const PressureConstraintSpec& spec,
        const EquationCompositionConfig& composition) {
    // The pressure predictor advances the same packed physical equations.
    // Its native WHICH method explicitly fuses these inputs without duplicating WHAT.
    Preset::installSingleFluid(system,SingleFluidPresetSpec{false},composition);
    system.addConstraint({"C_INCOMPRESSIBILITY","pressure/continuity constraint",
                          "pressureVelocityConsistency = 0","p"});
}

/// @brief Raw U/p composition for rho=rho0.  It intentionally does not
/// fabricate rhoE or a PerfectGas closure: execution is enabled only after a
/// dedicated constant-density predictor/corrector operator is bound.
///
/// STATE already supplies U/p; WHAT uses rho=rho0 to divide momentum and reduce continuity.
/// It never substitutes ddt(U) for variable-density ddt(rhoU).
/// Coupling preset 与 plan 由
/// pressure-coupling 模型贡献。
void addConstantDensityFluid(
        SystemCompositionBuilder& system,
        const EquationCompositionConfig& composition,
        bool includeDiffusion) {
    system.recordContribution("builtin.continuity","builtin incompressibility constraint");
    system.recordContribution("builtin.momentum","builtin Momentum equation");
    system.recordContribution("model.rhoConst","constant-density EOS closure");
    if (includeDiffusion) system.requireState("nu");
    auto momentumFormula=FormulaExpr::add(
        FormulaExpr::op("ddt",{FormulaExpr::symbol("U")},"momentum.time"),
        FormulaExpr::op("div",{FormulaExpr::symbol("momentumFlux")},
                        "momentum.convection"));
    if (includeDiffusion)
        momentumFormula=FormulaExpr::add(std::move(momentumFormula),
            FormulaExpr::op("diffusion",
                {FormulaExpr::symbol("nu"),FormulaExpr::symbol("U")},
                "momentum.diffusion"));
    system.addEquation({"momentum",std::move(momentumFormula),
        FormulaExpr::constantValue(0.0),{}});
    system.addEquation({"continuity",
        FormulaExpr::op("div",{FormulaExpr::symbol("U")},"continuity.constraint"),
        FormulaExpr::constantValue(0.0),{}});
    ExecutionScope call;
    call.kind=ExecutionKind::EquationCall;
    call.order=20;
    call.id="momentum.default";
    call.step={"momentum",{"U"}};
    call.origin={OriginKind::BuiltinPreset,"NavierStokes"};
    system.addExecution(std::move(call));
    system.addConstraint({"C_INCOMPRESSIBILITY","constant-density continuity",
                          "div(U) = 0","p"});
}

void addPhaseEquationPack(
        SystemCompositionBuilder& system,
        const std::string& phase,
        std::size_t phaseIndex, bool reference) {
    const std::string suffix = "."+phase;
    const std::string storage = "phase"+std::to_string(phaseIndex)+".";
    addState(system,"phaseMass"+suffix,"phase mass "+phase,1,
               StateRole::Transported,StorageBinding::ProviderDistributed,
               storage+"mass",0,phase);
    addState(system,"momentum"+suffix,"phase momentum "+phase,3,
               StateRole::Transported,StorageBinding::ProviderDistributed,
               storage+"momentum",0,phase);
    addState(system,"enthalpy"+suffix,"phase total enthalpy "+phase,1,
               StateRole::Transported,StorageBinding::ProviderDistributed,
               storage+"enthalpy",0,phase);
    const auto relations=eulerianPhaseRelations(phase,reference);
    for (const auto& equation:relations) system.addEquation(equation);
    for (const auto& item:std::vector<std::tuple<std::string,std::string,std::string,int>>{
            {"E_CONTINUITY"+suffix,"phaseMass"+suffix,"EulerianPhaseContinuity",10},
            {"momentum"+suffix,"momentum"+suffix,"EulerianPhaseMomentum",20},
            {"E_ENTHALPY"+suffix,"enthalpy"+suffix,"EulerianPhaseEnthalpy",60}}) {
        ExecutionScope call;call.kind=ExecutionKind::EquationCall;
        call.order=std::get<3>(item);call.step={std::get<0>(item),{std::get<1>(item)}};
        call.origin={OriginKind::Model,"Eulerian phase equations"};
        system.addExecution(std::move(call));
        system.bindNumerics({std::get<0>(item),std::get<2>(item)});
    }

}

void addSharedPressureConstraint(SystemCompositionBuilder& system) {
    system.addConstraint({
        "C_SHARED_PRESSURE","shared pressure relationship",
        "p.phase = p",""});
    system.addConstraint({
        "C_VOLUME_FRACTION","volume-fraction closure",
        "sum(alpha.phase) = 1",""});
    system.requestTransformation({
        "sharedPressureConstraint","PIMPLE shared-pressure transformation",
        110,true,{}});
}

void addEulerianEulerianTemplate(
        SystemCompositionBuilder& system,
        ResolvedSimulationSystem& resolved,
        const std::vector<std::string>& names,const std::string& referencePhase) {
    if (names.size() < 2) {
        throw std::runtime_error(
            "Resolved Eulerian-Eulerian system requires at least two phases.");
    }
    if (std::find(names.begin(),names.end(),referencePhase)==names.end())
        throw std::runtime_error("Eulerian reference phase must be explicitly declared in PhaseSystem.");
    system.recordContribution(
        "preset.eulerianEulerian","Eulerian-Eulerian equation preset");
    addState(system,"p","shared pressure",1,StateRole::Algebraic,
               StorageBinding::ProviderDistributed,"pressure",0,"pressure");
    for (std::size_t phase = 0; phase < names.size(); ++phase) {
        addPhaseEquationPack(system,names[phase],phase,names[phase]==referencePhase);
        for (const auto& item:std::vector<std::tuple<std::string,std::string,int>>{
                {"alpha","alpha",1},{"rho","density",1},{"U","velocity",3},
                {"h","primitiveEnthalpy",1},{"T","temperature",1}}) {
            addState(system,std::get<0>(item)+"."+names[phase],"recovered phase state",std::get<2>(item),
                StateRole::Derived,StorageBinding::ProviderDistributed,
                "phase"+std::to_string(phase)+"."+std::get<1>(item),0,names[phase]);
        }
        const std::string prefix = "phase"+std::to_string(phase)+".";
        resolved.runtime.workspaceRequirements.push_back({
            prefix+"momentumDiagonal",1,VariableLocation::EulerianCell,
            OwnershipKind::EulerianGlobalDof});
        resolved.runtime.workspaceRequirements.push_back({
            prefix+"volumeFaceFlux",3,VariableLocation::EulerianFace,
            OwnershipKind::CanonicalFace});
        resolved.runtime.workspaceRequirements.push_back({
            prefix+"massFaceFlux",3,VariableLocation::EulerianFace,
            OwnershipKind::CanonicalFace});
    }
    resolved.runtime.workspaceRequirements.push_back({
        "pressureCorrection",1,VariableLocation::EulerianCell,
        OwnershipKind::EulerianGlobalDof});
    system.addClosure("pressure reference row: pPrime(referenceCell)=0; homogeneous correction boundary closure follows existing pressure matrix");
    system.addClosure("referencePhaseRemainder = 1 - sum(non-reference alpha); existing recovery guards and EOS remain unchanged");
    system.addClosure("sumFaceAlphaSquaredOverDiagonal = sum_phase faceAverage(alpha^2 / momentumDiagonal)");
    system.addClosure("mixtureCompressibilityOverDt = sum(perfectGas alpha/p) / dt; pressureHistory = sum(perfectGas alpha/p)*(p-pOld)/dt");
    system.addClosure("effectivePhaseViscosity = alpha*(mu+mu_t); effectivePhaseConductivity = alpha*(k/Cp+eddyEnergyDiffusivity)");
    system.addClosure("materialPressureRate.phase = (p-pOld)/dt + U.phase dot gradient(p)");
    addSharedPressureConstraint(system);
}


} // namespace Compose
} // namespace SF::System
