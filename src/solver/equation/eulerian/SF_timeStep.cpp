/// @file SF_timeStep.cpp
/// @brief 双欧拉 CFL、扩散与相间刚性限制的局部时间步估计。

#include "solver/equation/eulerian/SF_equations.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace SF::EulerianEulerian {

double PhaseEquationAssembler::stableTimeStep(double cfl) const {
    if (!std::isfinite(cfl) || cfl <= 0.0) {
        throw std::runtime_error("Eulerian CFL must be finite and positive.");
    }
    const Field& field = system_.geometry();
    double limit = std::numeric_limits<double>::max();
    for (int cell : rowMap_.localCells) {
        int i = 0, j = 0, k = 0;
        field.getIJK(cell, i, j, k);
        for (const auto& state : system_.phases()) {
            for (int axis = 0; axis < 3; ++axis) {
                if (!Math::isDirectionActiveIndex(axis)) continue;
                int di = 0, dj = 0, dk = 0;
                Ops::offset(axis, 1, di, dj, dk);
                const double h = Ops::spacing(
                    field, i, j, k, i + di, j + dj, k + dk);
                const double nx =
                    (field.X(i+di,j+dj,k+dk)-field.X(i,j,k))/h;
                const double ny =
                    (field.Y(i+di,j+dj,k+dk)-field.Y(i,j,k))/h;
                const double nz =
                    (field.Z(i+di,j+dj,k+dk)-field.Z(i,j,k))/h;
                const double speed = std::abs(
                    state.primitive.velocity[0].values()[(size_t)cell]*nx
                    +state.primitive.velocity[1].values()[(size_t)cell]*ny
                    +state.primitive.velocity[2].values()[(size_t)cell]*nz);
                if (speed > 0.0) limit = std::min(limit,cfl*h/speed);
            }
        }
    }
    return limit;
}

double PhaseEquationAssembler::sourceTimeStep(double sourceCfl) const {
    if (!std::isfinite(sourceCfl)
        || sourceCfl <= 0.0 || sourceCfl > 1.0) {
        throw std::runtime_error(
            "Eulerian phaseSourceCFL must be in (0,1].");
    }
    double limit = std::numeric_limits<double>::max();
    for (int cell : rowMap_.localCells) {
        for (size_t phase = 0;
             phase < system_.phases().size(); ++phase) {
            const auto& state = system_.phases()[phase];
            const double mass =
                state.primary.phaseMass.values()[(size_t)cell];
            const double massRate = std::abs(
                system_.sources().mass[phase].values()[(size_t)cell]);
            if (massRate > 0.0) {
                limit = std::min(limit,sourceCfl*mass/massRate);
            }
            const double phaseEnthalpy =
                state.primary.phaseEnthalpy.values()[(size_t)cell];
            const double energyRate = std::abs(
                system_.sources().energy[phase].values()[(size_t)cell]);
            if (energyRate > 0.0) {
                limit = std::min(
                    limit,sourceCfl*phaseEnthalpy/energyRate);
            }
        }
        for (const auto& coupling :
             system_.sources().momentumCouplings) {
            const double drag =
                coupling.dragCoefficient.values()[(size_t)cell];
            if (drag <= 0.0) continue;
            const double first =
                system_.phases()[coupling.first].primary.phaseMass
                    .values()[(size_t)cell];
            const double second =
                system_.phases()[coupling.second].primary.phaseMass
                    .values()[(size_t)cell];
            limit = std::min(
                limit,sourceCfl*std::min(first,second)/drag);
        }
    }
    return limit;
}

void PhaseEquationAssembler::validateSourceStep(double dt) const {
    for (int cell : rowMap_.localCells) {
        for (size_t phase = 0;
             phase < system_.phases().size(); ++phase) {
            const auto& state = system_.phases()[phase];
            const double mass =
                state.primary.phaseMass.values()[(size_t)cell];
            const double massRate =
                system_.sources().mass[phase].values()[(size_t)cell];
            if (massRate < 0.0 && mass + dt*massRate <= 0.0) {
                throw std::runtime_error(
                    "Eulerian source would exhaust phaseMass at cell "
                    +std::to_string(cell)+", phase="+state.name
                    +"; reduce maxDeltaT.");
            }
            const double enthalpy =
                state.primary.phaseEnthalpy.values()[(size_t)cell];
            const double energyRate =
                system_.sources().energy[phase].values()[(size_t)cell];
            if (energyRate < 0.0
                && enthalpy + dt*energyRate <= 0.0) {
                throw std::runtime_error(
                    "Eulerian source would exhaust phaseEnthalpy at cell "
                    +std::to_string(cell)+", phase="+state.name
                    +"; reduce maxDeltaT.");
            }
        }
    }
}

} // namespace SF::EulerianEulerian
