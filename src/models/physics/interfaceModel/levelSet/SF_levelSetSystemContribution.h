#pragma once

/// @file SF_levelSetSystemContribution.h
/// @brief Level-set-owned mathematical system contribution.

#include "core/system/SF_systemContribution.h"

namespace SF::Physics::InterfaceModels::LevelSetContribution {

struct Spec {
    bool ghostFluid = false;
    int advectionOrder = 5;
    int reinitializationOrder = 5;
    int reinitializationSteps = 0;
    double pseudoTimeStep = 0.0;
    double wenoEpsilon = 1e-6;
    double wenoPower = 2.0;
    double signSmoothingFactor = 1.0;
    bool continuousSurfaceForce = false;
    double surfaceTension = 0.0;
    double interfaceThickness = 0.0;
};

void contribute(System::SystemContribution& system, const Spec& spec);

} // namespace SF::Physics::InterfaceModels::LevelSetContribution
