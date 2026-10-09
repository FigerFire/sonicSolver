#pragma once
#include "SF_phaseProperties.h"
#include "core/system/SF_systemContribution.h"
namespace SF::Physics::Multiphase {
/// The physical model expands native mathematics and storage; no timestep ownership.
void contributeHomogeneous(System::SystemContribution& system,const MultiPhaseConfig& config);
void contributeMixture(System::SystemContribution& system,const MultiPhaseConfig& config);
}
