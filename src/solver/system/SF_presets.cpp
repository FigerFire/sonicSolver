#include "SF_builtinState.h"
/// @file SF_presets.cpp
/// @brief COMPOSE — built-in 物理方程 bundle 的四模块贡献。
///
/// 默认 equation calls 与 numerical bindings 随数学贡献声明；全局 coupling/time
/// 仍来自独立显式请求，Equation pack 不选择 loop topology。

#include "SF_presets.h"
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
ExecutionProgram explicitProgram(std::vector<ExecutionScope> calls) {
    ExecutionProgram program;program.recipeScopes=true;
    program.root.children=std::move(calls);
    ExecutionScope commit;commit.kind=ExecutionKind::Commit;
    commit.id="physicalStep.commit";commit.order=1000000;
    commit.origin={OriginKind::BuiltinPreset,"builtin explicit"};
    program.root.children.push_back(std::move(commit));
    return program;
}

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
    if(binding==StorageBinding::TransientWorkspace)
        unknown.initializationRequired=unknown.boundaryRequired=unknown.outputEligible=unknown.runtimeStorageRequired=unknown.restartEligible=false;
    unknown.signature=EulerianMathematics::signature(unknown.id);
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
                          "pressureVelocityConsistency = 0","p",{},{"pSimple","correctU","correctFluxp"}});
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
        bool includeDiffusion,double dynamicViscosity) {
    system.recordContribution("builtin.continuity","builtin incompressibility constraint");
    system.recordContribution("model.rhoConst","constant-density EOS closure");
    using E=FormulaExpr;
    const auto u=BuiltinStateCatalog().expression("U");
    const auto p=BuiltinStateCatalog().expression("p");
    const auto rho=BuiltinStateCatalog().expression("rho");
    const auto& predictor=composition.predictorEquation;
    if(predictor.empty()) {
        system.recordContribution("builtin.PressureMomentum","typed constant-density momentum declaration");
        auto spatial=E::negate(E::op(Math::Op::Div,{E::op(Math::Op::Outer,{u,u})},"momentum.convection"));
        if(includeDiffusion) {
            if(composition.thermoDynamics.transport!="const")throw std::runtime_error("Cell explicit constant-density preset needs a constant transport coefficient; no variable-property response is implemented");
            auto nu=BuiltinStateCatalog().at("nu");nu.producer.reset();nu.evaluation=StateEvaluation::Direct;
            nu.constantValue=dynamicViscosity/composition.thermoDynamics.constantDensity;
            nu.storageBinding=StorageBinding::SpecializedExecutor;nu.storageKey="nuConst";
            nu.runtimeStorageRequired=nu.initializationRequired=nu.boundaryRequired=nu.outputEligible=nu.restartEligible=false;
            system.addState(nu);
            spatial=E::add(std::move(spatial),E::op(Math::Op::Diffusion,{nu.expression(),u},pressureMomentumDiffusion));
        }
        const auto force=E::negate(E::divide(E::op(Math::Op::Grad,{p},"momentum.pressureGradient"),rho));
        system.addEquation({pressureMomentumEquation,E::op(Math::Op::Ddt,{u},"momentum.time"),E::add(std::move(spatial),force),{}});
    }
    system.addEquation({"continuity",
        FormulaExpr::op(Math::Op::Div,{BuiltinStateCatalog().expression("U")},"continuity.constraint"),
        FormulaExpr::constantValue(0.0),{}});
    ExecutionScope call;
    call.kind=ExecutionKind::EquationCall;
    call.order=20;
    call.id="momentum.default";
    call.step={predictor.empty()?pressureMomentumEquation:predictor,{"U"}};
    call.origin={OriginKind::BuiltinPreset,"NavierStokes"};
    system.pressureInputs=pressureCouplingInputs(call.step,PressurePredictorKind::Velocity);
    system.addExecution(std::move(call));
    system.addConstraint({"C_INCOMPRESSIBILITY","constant-density continuity",
                          "div(U) = 0","p",{},{"continuity","pSimple","correctU","correctFluxp"}});
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
        "p.phase = p","",{},{"E_SHARED_PRESSURE","correctSharedP"}});
    std::vector<std::string> phaseBalances;
    for (const auto& equation:system.rawSystem().registry.entries())
        if (equation.id.rfind("E_CONTINUITY.",0)==0) phaseBalances.push_back(equation.id);
    system.addConstraint({"C_VOLUME_FRACTION","volume-fraction closure",
        "sum(alpha.phase) = 1","",{},std::move(phaseBalances)});
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
    addState(system,"deltaP","shared pressure increment",1,StateRole::Algebraic,StorageBinding::TransientWorkspace,"pressureCorrection");
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
    std::vector<std::string> phaseBalances,phaseStressAndHeat,phaseHeat;
    for (const auto& name:names) {
        phaseBalances.push_back("E_CONTINUITY."+name);
        phaseStressAndHeat.push_back("momentum."+name);
        phaseStressAndHeat.push_back("E_ENTHALPY."+name);
        phaseHeat.push_back("E_ENTHALPY."+name);
    }
    system.addClosure("pressure reference row: deltaP(referenceCell)=0; homogeneous correction boundary closure follows existing pressure matrix",{"E_SHARED_PRESSURE"});
    system.addClosure("referencePhaseRemainder = 1 - sum(non-reference alpha); existing recovery guards and EOS remain unchanged",phaseBalances);
    system.addClosure("sumFaceAlphaSquaredOverDiagonal = sum_phase faceAverage(alpha^2 / momentumDiagonal)",{"E_SHARED_PRESSURE"});
    system.addClosure("mixtureCompressibilityOverDt = sum(perfectGas alpha/p) / dt; pressureHistory = sum(perfectGas alpha/p)*(p-pOld)/dt",{"E_SHARED_PRESSURE"});
    system.addClosure("effectivePhaseViscosity = alpha*(mu+mu_t); effectivePhaseConductivity = alpha*(k/Cp+eddyEnergyDiffusivity)",phaseStressAndHeat);
    system.addClosure("materialPressureRate.phase = (p-pOld)/dt + U.phase dot gradient(p)",phaseHeat);
    addSharedPressureConstraint(system);
}


} // namespace Compose
} // namespace SF::System
