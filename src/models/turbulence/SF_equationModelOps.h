#pragma once

/// @file SF_equationModelOps.h
/// @brief Eulerian 湍流方程闭式共享的几何与微分算子。

#include "SF_equationSystem.h"

#include <array>
#include <string>

namespace SF::Turbulence::ModelOps {

bool isPhysical(const Field& geometry, int i, int j, int k);

double strainSquared(
    const Field& geometry,
    const Physics::PhaseSystems::PhaseState& phase,
    int i, int j, int k);

double vorticityMagnitude(
    const Field& geometry,
    const Physics::PhaseSystems::PhaseState& phase,
    int i, int j, int k);

double filterWidth(
    const Field& geometry,
    int i, int j, int k,
    double scale);

void clearEquationCell(
    TransportEquationState& equation,
    int cell);

void requirePositiveState(
    const std::string& model,
    const std::string& phase,
    const std::string& secondName,
    double density,
    double phaseMass,
    double kineticEnergy,
    double secondVariable,
    double kineticEnergyFloor,
    double secondFloor,
    int i, int j, int k);

void buildWallDistance(
    const Physics::PhaseSystems::PhaseSystem& system,
    const FDM::TurbulenceConfig& config,
    PhaseEquationState& state);

} // namespace SF::Turbulence::ModelOps
