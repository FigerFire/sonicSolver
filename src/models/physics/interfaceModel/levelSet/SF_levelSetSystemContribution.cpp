/// @file SF_levelSetSystemContribution.cpp
/// @brief Level-set unknown, equation, and closure contribution.

#include "SF_levelSetSystemContribution.h"

#include "core/system/SF_systemContribution.h"

#include <utility>

namespace SF::Physics::InterfaceModels::LevelSetContribution {

void contribute(
        System::SystemContribution& system, const Spec& spec) {
    system.recordContribution("model.levelSet","level-set equations");
    System::UnknownDescriptor phi;
    phi.id = "phi";
    phi.name = "level-set geometry";
    phi.components = 1;
    phi.shape = System::ValueShape::Scalar;
    phi.role = System::UnknownRole::Transported;
    phi.storageBinding = System::StorageBinding::NamedDistributed;
    phi.storageKey = "phi";
    phi.nameSpace = "interface";
    system.addUnknown(std::move(phi));
    system.addEquation(
        {"E_LEVEL_SET","level-set advection/reinitialization",
         "geometry transport",{"phi"}},
        Equation::named("E_LEVEL_SET",
            Equation::ddt({"phi"}) + Equation::div({"levelSetFlux"})
                == Equation::Symbol{"reinitialization"}));
    system.addClosure("surface normal and curvature from phi");
    if (spec.ghostFluid) system.addClosure("ghost-fluid interface closure");
}

} // namespace SF::Physics::InterfaceModels::LevelSetContribution
