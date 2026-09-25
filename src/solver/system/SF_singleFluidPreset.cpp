/// @file SF_singleFluidPreset.cpp
/// @brief Authoritative density single-fluid Navier-Stokes composition.

#include "SF_singleFluidPreset.h"

#include <utility>

namespace SF::System::Preset {
namespace {

UnknownDescriptor packedUnknown(
        std::string id, std::string name, int components, int offset) {
    UnknownDescriptor unknown;
    unknown.id = std::move(id);
    unknown.name = std::move(name);
    unknown.components = components;
    unknown.shape = components == 1 ? ValueShape::Scalar : ValueShape::Vector;
    unknown.role = UnknownRole::Primary;
    unknown.storageBinding = StorageBinding::PackedDistributed;
    unknown.storageKey = "conservative";
    unknown.componentOffset = offset;
    unknown.nameSpace = "fluid";
    return unknown;
}

} // namespace

void installSingleFluid(
        SystemCompositionBuilder& system,
        const SingleFluidPresetSpec& spec) {
    system.recordContribution(
        "builtin.singleFluidNavierStokes",
        "density single-fluid Navier-Stokes equation preset");
    system.addUnknown(packedUnknown("rho","density",1,0));
    system.addUnknown(packedUnknown("rhoU","momentum",3,1));
    system.addUnknown(packedUnknown("rhoE","total energy",1,4));

    system.addEquation({"E_MASS","mass","conservation",{"rho"}},
        Equation::named("E_MASS",
            Equation::ddt({"rho"}) + Equation::div({"massFlux"})
                == Equation::Symbol{"zero"}));

    auto momentum = Equation::ddt({"rhoU"})
        + Equation::div({"momentumFlux"});
    auto energy = Equation::ddt({"rhoE"})
        + Equation::div({"energyFlux"});
    if (spec.diffusion) {
        momentum = std::move(momentum)
            + Equation::diffusion({"mu"},{"U"});
        energy = std::move(energy)
            + Equation::diffusion({"conductivity"},{"T"});
    }

    system.addEquation({"E_MOMENTUM","momentum","conservation",{"rhoU"}},
        Equation::named("E_MOMENTUM",
            std::move(momentum) == Equation::Symbol{"zero"}));
    system.addEquation({"E_ENERGY","total energy","conservation",{"rhoE"}},
        Equation::named("E_ENERGY",
            std::move(energy) == Equation::Symbol{"zero"}));
}

} // namespace SF::System::Preset
