/// @file SF_turbulenceSystemContribution.cpp
/// @brief Turbulence unknown, equation, closure, and schedule contributions.

#include "SF_turbulenceSystemContribution.h"
#include "core/system/SF_scheduleIds.h"

#include "core/system/SF_systemContribution.h"

#include <stdexcept>
#include <utility>

namespace SF::Turbulence {
namespace {

void addScalar(
        System::SystemContribution& system,
        const std::string& id,
        const std::string& name,
        const std::string& storage,
        const std::string& nameSpace) {
    System::UnknownDescriptor unknown;
    unknown.id = id;
    unknown.name = name;
    unknown.components = 1;
    unknown.shape = System::ValueShape::Scalar;
    unknown.role = System::UnknownRole::Transported;
    unknown.storageBinding = System::StorageBinding::NamedDistributed;
    unknown.storageKey = storage;
    unknown.nameSpace = nameSpace;
    system.addUnknown(std::move(unknown));
}

} // namespace

void contribute(
        System::SystemContribution& system,
        const SystemContributionSpec& spec) {
    system.recordContribution(
        "model.turbulence","turbulence equation/closure contribution");
    const bool kEpsilon = spec.model == "kEpsilon";
    const bool kOmega = spec.model == "kOmegaSST";
    if (!kEpsilon && !kOmega) {
        system.addClosure(
            spec.model.empty() ? "turbulence closure" : "turbulence: "+spec.model);
        return;
    }

    std::vector<std::string> phases = spec.phases;
    if (!spec.eulerian) phases = {""};
    else if (phases.empty()) {
        throw std::runtime_error(
            "Eulerian transported turbulence equations require explicit phases.");
    }
    for (std::size_t phaseIndex = 0; phaseIndex < phases.size(); ++phaseIndex) {
        const std::string& phase = phases[phaseIndex];
        const std::string suffix = phase.empty() ? "" : "."+phase;
        const std::string k = "k"+suffix;
        const std::string second = (kEpsilon ? "epsilon" : "omega")+suffix;
        const std::string kEquation = "E_TURB_K"+suffix;
        const std::string secondEquation =
            (kEpsilon ? "E_TURB_EPSILON" : "E_TURB_OMEGA")+suffix;
        const std::string storage = phase.empty()
            ? "" : "phase"+std::to_string(phaseIndex)+".";
        addScalar(system,k,"turbulent kinetic energy"+
            (phase.empty() ? "" : " "+phase),storage+"k",
            phase.empty() ? "turbulence" : phase);
        addScalar(system,second,
            std::string(kEpsilon ? "turbulence dissipation "
                                 : "specific dissipation ")+phase,
            storage+(kEpsilon ? "epsilon" : "omega"),
            phase.empty() ? "turbulence" : phase);
        system.addEquation(
            {kEquation,"turbulent kinetic energy transport","conservation",{k}},
            Equation::named(kEquation,
                Equation::ddt({k}) + Equation::div({"flux."+k})
                    + Equation::diffusion({"diffusivity."+k},{k})
                    == Equation::Symbol{"source."+k}));
        system.addEquation(
            {secondEquation,kEpsilon ? "epsilon transport" : "omega transport",
             "conservation",{second}},
            Equation::named(secondEquation,
                Equation::ddt({second}) + Equation::div({"flux."+second})
                    + Equation::diffusion(
                        {"diffusivity."+second},{second})
                    == Equation::Symbol{"source."+second}));
        if (spec.eulerian) {
            system.extendExecutionPolicy(System::kSharedPressureScheduleId,kEquation,k);
            system.extendExecutionPolicy(System::kSharedPressureScheduleId,secondEquation,second);
        }
    }
    system.addClosure("mu_t from "+spec.model);
}

} // namespace SF::Turbulence
