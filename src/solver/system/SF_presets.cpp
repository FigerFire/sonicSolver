/// @file SF_presets.cpp
/// @brief COMPOSE — built-in 方程族 preset 的 composition（WHAT only）。
///
/// 这里只声明未知量、方程与约束。Coupling preset / plan / backend / time
/// recipe 不由方程 preset 选择。

#include "SF_systemBuilder.h"
#include "SF_equationContribution.h"
#include "SF_singleFluidPreset.h"

#include "SF_config.h"

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

namespace SF::System {
namespace Compose {
namespace {

bool usesEquation(
        const EquationCompositionConfig& composition,
        const std::string& name) {
    return std::find(composition.equations.begin(),composition.equations.end(),
                     name) != composition.equations.end();
}

void addUnknown(
        SystemCompositionBuilder& system,
        std::string id,
        std::string name,
        int components = 1,
        UnknownRole role = UnknownRole::Primary,
        StorageBinding binding = StorageBinding::SpecializedExecutor,
        std::string storageKey = {},
        int componentOffset = 0,
        std::string nameSpace = {}) {
    UnknownDescriptor unknown;
    unknown.id = std::move(id);
    unknown.name = std::move(name);
    unknown.components = components;
    unknown.shape = components == 1 ? ValueShape::Scalar : ValueShape::Vector;
    unknown.role = role;
    unknown.storageBinding = binding;
    unknown.storageKey = std::move(storageKey);
    unknown.componentOffset = componentOffset;
    unknown.nameSpace = std::move(nameSpace);
    system.addUnknown(std::move(unknown));
}

void addEquation(
        SystemCompositionBuilder& system,
        EquationDescriptor descriptor,
        Equation::Definition definition) {
    system.addEquation(std::move(descriptor),std::move(definition));
}

} // namespace

void addPressureConstraintFluid(
        SystemCompositionBuilder& system,
        const PressureConstraintSpec& spec) {
    system.recordContribution("builtin.navierStokes","core Navier-Stokes preset");
    addUnknown(system,"rho","density",1,UnknownRole::Primary,
               StorageBinding::PackedDistributed,"conservative",0,"fluid");
    addUnknown(system,"U","velocity",3,UnknownRole::Primary,
               StorageBinding::PackedDistributed,"conservative",1,"fluid");
    addUnknown(system,"rhoE","total energy",1,UnknownRole::Primary,
               StorageBinding::PackedDistributed,"conservative",4,"fluid");
    addUnknown(system,"p","thermodynamic pressure",1,
               UnknownRole::Derived,StorageBinding::SpecializedExecutor,
               "thermodynamicClosure",0,"pressure");
    addEquation(system,{"E_MASS","mass","conservation",{"rho"}},
        Equation::named("E_MASS",
            Equation::ddt({"rho"}) + Equation::div({"massFlux"})
                == Equation::Symbol{"zero"}));
    addEquation(system,
        {"E_MOMENTUM","momentum predictor","conservation",{"U"}},
        Equation::named("E_MOMENTUM",
            Equation::ddt({"U"}) + Equation::div({"momentumFlux"})
                == Equation::Symbol{"zero"}));
    addEquation(system,
        {"E_ENERGY","total energy","conservation",{"rhoE"}},
        Equation::named("E_ENERGY",
            Equation::ddt({"rhoE"}) + Equation::div({"energyFlux"})
                == Equation::Symbol{"zero"}));
    ConstraintDescriptor pressureConstraint{
        "C_INCOMPRESSIBILITY","pressure/continuity constraint",
        "pressureVelocityConsistency = 0","p"};
    system.addConstraint(std::move(pressureConstraint));
}

/// @brief Raw U/p composition for rho=rho0.  It intentionally does not
/// fabricate rhoE or a PerfectGas closure: execution is enabled only after a
/// dedicated constant-density predictor/corrector operator is bound.
///
/// 它只声明 WHAT（守恒量 + div 约束 + 乘子）；coupling preset 与 plan 由
/// pressure-coupling 模型贡献。
void addConstantDensityFluid(
        SystemCompositionBuilder& system,
        const EquationCompositionConfig& composition,
        bool includeDiffusion) {
    system.recordContribution("builtin.continuity","builtin Continuity equation");
    system.recordContribution("builtin.momentum","builtin Momentum equation");
    system.recordContribution("model.rhoConst","constant-density EOS closure");
    addUnknown(system,"rho","constant density closure",1,UnknownRole::Derived,
               StorageBinding::SpecializedExecutor,"rhoConst",0,"fluid");
    addUnknown(system,"U","velocity",3,UnknownRole::Primary,
               StorageBinding::PackedDistributed,"conservative",1,"fluid");
    addUnknown(system,"p","pressure multiplier",1,UnknownRole::Multiplier,
               StorageBinding::NamedDistributed,"pressure",0,"pressure");
    addEquation(system,{"E_CONTINUITY","continuity","constraint",{"U"}},
        Equation::named("E_CONTINUITY",
            Equation::div({"U"}) == Equation::Symbol{"zero"}));
    auto momentum = Equation::ddt({"U"})
        + Equation::div({"momentumFlux"});
    if (includeDiffusion) {
        momentum = std::move(momentum) + Equation::diffusion({"nu"},{"U"});
    }
    addEquation(system,{"E_MOMENTUM","momentum","conservation",{"U"}},
        Equation::named("E_MOMENTUM",
            std::move(momentum) == Equation::Symbol{"zero"}));
    system.addConstraint({"C_INCOMPRESSIBILITY","constant-density continuity",
                          "div(U) = 0","p"});
}

void addPhaseEquationPack(
        SystemCompositionBuilder& system,
        const std::string& phase,
        std::size_t phaseIndex) {
    const std::string suffix = "."+phase;
    const std::string storage = "phase"+std::to_string(phaseIndex)+".";
    addUnknown(system,"phaseMass"+suffix,"phase mass "+phase,1,
               UnknownRole::Transported,StorageBinding::NamedDistributed,
               storage+"mass",0,phase);
    addUnknown(system,"momentum"+suffix,"phase momentum "+phase,3,
               UnknownRole::Transported,StorageBinding::NamedDistributed,
               storage+"momentum",0,phase);
    addUnknown(system,"enthalpy"+suffix,"phase total enthalpy "+phase,1,
               UnknownRole::Transported,StorageBinding::NamedDistributed,
               storage+"enthalpy",0,phase);
    addEquation(system,{
        "E_CONTINUITY"+suffix,"phase continuity "+phase,
        "conservation",{"phaseMass"+suffix}},
        Equation::named("E_CONTINUITY"+suffix,
            Equation::ddt({"phaseMass"+suffix})
                + Equation::div({"phaseMassFlux"+suffix})
                == Equation::Symbol{"phaseMassSources"+suffix}));
    addEquation(system,{
        "E_MOMENTUM"+suffix,"phase momentum "+phase,
        "conservation",{"momentum"+suffix}},
        Equation::named("E_MOMENTUM"+suffix,
            Equation::ddt({"momentum"+suffix})
                + Equation::div({"phaseMomentumFlux"+suffix})
                + Equation::gradient({"alphaPressure"+suffix})
                + Equation::diffusion(
                    {"effectivePhaseViscosity"+suffix},{"U"+suffix})
                == Equation::Symbol{"phaseMomentumSources"+suffix}));
    addEquation(system,{
        "E_ENTHALPY"+suffix,"phase enthalpy "+phase,
        "conservation",{"enthalpy"+suffix}},
        Equation::named("E_ENTHALPY"+suffix,
            Equation::ddt({"enthalpy"+suffix})
                + Equation::div({"phaseEnthalpyFlux"+suffix})
                + Equation::diffusion(
                    {"effectivePhaseConductivity"+suffix},{"h"+suffix})
                == Equation::Symbol{"phaseEnthalpySources"+suffix}));
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
        const std::vector<std::string>& names) {
    if (names.size() < 2) {
        throw std::runtime_error(
            "Resolved Eulerian-Eulerian system requires at least two phases.");
    }
    system.recordContribution(
        "preset.eulerianEulerian","Eulerian-Eulerian equation preset");
    system.recordContribution(
        "preset.coupling.PIMPLE","PIMPLE shared-pressure preset");
    addUnknown(system,"p","shared pressure",1,UnknownRole::Algebraic,
               StorageBinding::NamedDistributed,"pressure",0,"pressure");
    for (std::size_t phase = 0; phase < names.size(); ++phase) {
        addPhaseEquationPack(system,names[phase],phase);
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
    addSharedPressureConstraint(system);
}


} // namespace Compose
} // namespace SF::System
