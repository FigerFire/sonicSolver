/// @file SF_singleFluidPreset.cpp
/// @brief Authoritative density single-fluid Navier-Stokes composition.

#include "SF_singleFluidPreset.h"

#include <algorithm>
#include <utility>
#include <tuple>

namespace SF::System::Preset {
namespace {

FormulaExpr transient(const char* symbol,const char* occurrence) {
    return FormulaExpr::op("ddt",{FormulaExpr::symbol(symbol)},occurrence);
}
FormulaExpr divergence(const char* flux,const char* occurrence) {
    return FormulaExpr::op("div",{FormulaExpr::symbol(flux)},occurrence);
}
FormulaExpr diffusion(const char* coefficient,const char* unknown,
                      const char* occurrence) {
    return FormulaExpr::op("diffusion",
        {FormulaExpr::symbol(coefficient),FormulaExpr::symbol(unknown)},
        occurrence);
}

} // namespace

void installSingleFluid(
        SystemCompositionBuilder& system,
        const SingleFluidPresetSpec& spec,
        const EquationCompositionConfig& composition) {
    system.recordContribution(
        "builtin.singleFluidNavierStokes",
        "single-fluid conservation balance equations");
    for (const auto* symbol:{"U","p","T"}) system.requireState(symbol);
    if (spec.diffusion) {
        system.requireState("mu");
        system.requireState("conductivity");
    }

    const auto selected=[&](const char* name) {
        return composition.equations.empty() || std::find(composition.equations.begin(),composition.equations.end(),name)!=composition.equations.end();
    };
    if (selected("Continuity")) system.addEquation({"continuity",FormulaExpr::add(transient("rho","mass.time"),
        divergence("massFlux","mass.convection")),FormulaExpr::constantValue(0.0),{}});
    auto momentum=FormulaExpr::add(transient("rhoU","momentum.time"),
        divergence("momentumFlux","momentum.convection"));
    auto energy=FormulaExpr::add(transient("rhoE","energy.time"),
        divergence("energyFlux","energy.convection"));
    if (spec.diffusion) {
        momentum=FormulaExpr::add(std::move(momentum),diffusion("mu","U","momentum.diffusion"));
        energy=FormulaExpr::add(std::move(energy),diffusion("conductivity","T","energy.diffusion"));
    }
    if (selected("Momentum")) system.addEquation({"momentum",std::move(momentum),FormulaExpr::constantValue(0.0),{}});
    if (selected("Energy")) system.addEquation({"energy",std::move(energy),FormulaExpr::constantValue(0.0),{}});
    for (const auto& item:std::vector<std::tuple<const char*,const char*,int>>{
            {"continuity","rho",10},{"momentum","rhoU",20},{"energy","rhoE",60}}) {
        if (!system.rawSystem().registry.contains(std::get<0>(item))) continue;
        const auto target=std::string(std::get<0>(item))=="momentum"
            && system.rawSystem().state.isSolution("U") ? "U" : std::get<1>(item);
        if (!system.rawSystem().state.isSolution(target))
            throw std::runtime_error(std::string("No selected solution STATE target for equation ")+std::get<0>(item));
        ExecutionScope call;
        call.kind=ExecutionKind::EquationCall;
        call.order=std::get<2>(item);
        call.step={std::get<0>(item),{target}};
        call.origin={OriginKind::BuiltinPreset,"NavierStokes"};
        system.addExecution(std::move(call));
        system.bindNumerics({std::get<0>(item),"ConservativeResidual"});
    }
    ExecutionScope commit;
    commit.kind=ExecutionKind::Commit;
    commit.id="physicalStep.commit";
    commit.order=1000000;
    commit.origin={OriginKind::BuiltinPreset,"NavierStokes"};
    system.addExecution(std::move(commit));
}

} // namespace SF::System::Preset
