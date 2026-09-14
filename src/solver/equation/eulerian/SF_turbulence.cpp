/// @file SF_turbulence.cpp
/// @brief 将各相湍流输运方程自由耦合到双欧拉外迭代。

#include "solver/equation/eulerian/SF_equations.h"

#include "solver/equation/eulerian/SF_transport.h"

#include <algorithm>
#include <stdexcept>

namespace SF::EulerianEulerian {

void PhaseEquationAssembler::attachTurbulence(
        Turbulence::EquationSystem* turbulence) {
    turbulence_ = turbulence;
    turbulenceSolvers_.clear();
    if (!turbulence_ || !turbulence_->hasTransportEquations()) return;
    for (auto& state : turbulence_->states()) {
        for (const auto* equation :
             turbulence_->transportEquations(state)) {
            (void)equation;
            turbulenceSolvers_.push_back(
                std::make_unique<LinearAlgebra::SolverSession>(
                    config_.turbulence));
        }
    }
}

void PhaseEquationAssembler::resetTurbulenceDiagnostics() {
    turbulenceIterations_ = 0;
    turbulenceResidual_ = 0.0;
}

LinearAlgebra::ReuseStatistics
PhaseEquationAssembler::turbulenceReuseStatistics() const {
    LinearAlgebra::ReuseStatistics total;
    for (const auto& solver : turbulenceSolvers_) {
        const auto& value = solver->statistics();
        total.structureRebuilds += value.structureRebuilds;
        total.coefficientUpdates += value.coefficientUpdates;
        total.rhsUpdates += value.rhsUpdates;
        total.preconditionerRefreshes +=
            value.preconditionerRefreshes;
        total.solves += value.solves;
    }
    return total;
}

void PhaseEquationAssembler::solveTurbulence(double timeStep) {
    if (!turbulence_ || !turbulence_->hasTransportEquations()) return;
    const Field& field = system_.geometry();
    ScalarField diagonal;
    diagonal.setupLike(field, "turbulenceDiagonal");
    size_t solverIndex = 0;

    for (auto& phaseState : turbulence_->states()) {
        const size_t phase = phaseState.phaseIndex;
        auto solveEquation = [&](
                Turbulence::TransportEquationState& equation) {
            auto matrix = assembleTransportEquation(
                field, rowMap_, equation.variable,
                equation.previousConserved,
                system_.phases()[phase].primary.phaseMass,
                equation.diffusivity, equation.explicitSource,
                workspace_.faceFlux[phase].mass, timeStep,
                system_.config().eulerianEulerian.axisymmetric,
                system_.config().eulerianEulerian.radialCoordinate,
                diagonal, &equation.implicitSink);
            LinearAlgebra::SolveResult result;
            try {
                result =
                    turbulenceSolvers_.at(solverIndex)->solve(matrix);
            } catch (const std::exception& error) {
                throw std::runtime_error(
                    "Eulerian turbulence solve failed for '"
                    + equation.variable.name() + "': " + error.what());
            }
            ++solverIndex;
            turbulenceIterations_ += result.iterations;
            turbulenceResidual_ = std::max(
                turbulenceResidual_, result.relativeResidual);
            for (size_t row = 0; row < rowMap_.localCells.size(); ++row) {
                equation.variable.values()
                    [(size_t)rowMap_.localCells[row]] =
                        result.solution[row];
            }
        };
        for (auto* equation :
             turbulence_->transportEquations(phaseState)) {
            solveEquation(*equation);
        }
    }

    turbulence_->applyBoundary(system_);
    turbulence_->validate(system_, "linear solve");
    turbulence_->prepare(system_);
}

} // namespace SF::EulerianEulerian
