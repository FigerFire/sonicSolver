/// @file SF_phase.cpp
/// @brief 双欧拉每相连续、动量和能量离散系统的组装求解。

#include "solver/equation/eulerian/SF_equations.h"

#include "SF_linearSystem.h"
#include "solver/equation/eulerian/SF_transport.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace SF::EulerianEulerian {
namespace {

double faceDivergence(const Field& field,
                      const std::vector<double>& flux,
                      int i, int j, int k,
                      bool axisymmetric, int radialCoordinate) {
    double divergence = 0.0;
    for (int axis = 0; axis < 3; ++axis) {
        if (!Math::isDirectionActiveIndex(axis)) continue;
        int di = 0, dj = 0, dk = 0;
        Ops::offset(axis, -1, di, dj, dk);
        const size_t high = PhaseFaceFlux::index(
            field, axis, i, j, k);
        const size_t low = PhaseFaceFlux::index(
            field, axis, i + di, j + dj, k + dk);
        divergence += flux[high] - flux[low];
    }
    return Ops::inverseCellVolume(
        field, i, j, k, axisymmetric, radialCoordinate) * divergence;
}

LinearAlgebra::DistributedRowMap makePhaseRowMap(
        Physics::PhaseSystems::PhaseSystem& system,
        FDM::IExecutionRuntime* runtime) {
    const bool axisymmetric =
        system.config().eulerianEulerian.axisymmetric;
    const int radialCoordinate =
        system.config().eulerianEulerian.radialCoordinate;
    return LinearAlgebra::makeDistributedRowMap(
        system.geometry(), runtime,
        [&system, axisymmetric, radialCoordinate](int i, int j, int k) {
            return Ops::solved(
                system.geometry(), i, j, k,
                axisymmetric, radialCoordinate);
        });
}

} // namespace

void PhaseEquationAssembler::requireTerm(
        const Equation::AssemblyPlan& plan,
        Equation::TermKind kind,
        const char* operation) const {
    if (!plan.contains(kind)) {
        throw std::runtime_error(
            "Eulerian runtime operation '"+std::string(operation)
            +"' requested an operator absent from resolved equation '"
            +plan.definition->name+"'.");
    }
}

void PhaseEquationAssembler::bindAssemblyPlans(
        const Equation::AssemblyPlanRegistry& plans) {
    if (!phasePlans_.empty() || pressurePlan_) {
        throw std::runtime_error(
            "Eulerian AssemblyPlans are already bound.");
    }
    phasePlans_.reserve(system_.phases().size());
    for (const auto& phase : system_.phases()) {
        const std::string suffix = "."+phase.name;
        PhasePlans bound{
            &plans.at("E_CONTINUITY"+suffix),
            &plans.at("E_MOMENTUM"+suffix),
            &plans.at("E_ENTHALPY"+suffix)};
        requireTerm(*bound.continuity,Equation::TermKind::Transient,
                    "continuity time derivative");
        requireTerm(*bound.continuity,Equation::TermKind::Divergence,
                    "continuity face divergence");
        requireTerm(*bound.momentum,Equation::TermKind::Transient,
                    "momentum time derivative");
        requireTerm(*bound.momentum,Equation::TermKind::Divergence,
                    "momentum convection");
        requireTerm(*bound.momentum,Equation::TermKind::Gradient,
                    "momentum pressure gradient");
        requireTerm(*bound.momentum,Equation::TermKind::Diffusion,
                    "momentum diffusion");
        requireTerm(*bound.enthalpy,Equation::TermKind::Transient,
                    "enthalpy time derivative");
        requireTerm(*bound.enthalpy,Equation::TermKind::Divergence,
                    "enthalpy convection");
        requireTerm(*bound.enthalpy,Equation::TermKind::Diffusion,
                    "enthalpy diffusion");
        phasePlans_.push_back(bound);
    }
    pressurePlan_ = &plans.at("E_SHARED_PRESSURE");
    requireTerm(*pressurePlan_,Equation::TermKind::Constraint,
                "shared-pressure correction");
}

PhaseEquationAssembler::PhaseEquationAssembler(
        Physics::PhaseSystems::PhaseSystem& system,
        const FDM::PressureCorrectionConfig& config,
        PhaseSolverWorkspace& workspace)
    : system_(system), config_(config), workspace_(workspace),
      rowMap_(makePhaseRowMap(system, nullptr)),
      pressureSolver_(config.linear.pressure) {
    const size_t phaseCount = system_.phases().size();
    momentumSolvers_.reserve(3 * phaseCount);
    energySolvers_.reserve(phaseCount);
    for (size_t phase = 0; phase < phaseCount; ++phase) {
        for (int component = 0; component < 3; ++component) {
            momentumSolvers_.push_back(
                std::make_unique<LinearAlgebra::SolverSession>(
                    config_.linear.momentum));
        }
        energySolvers_.push_back(
            std::make_unique<LinearAlgebra::SolverSession>(
                config_.linear.energy));
    }
}

void PhaseEquationAssembler::refreshRowMap() {
    rowMap_ = makePhaseRowMap(system_, runtime_);
}

void PhaseEquationAssembler::initializeMomentumDiagonal(double dt) {
    if (phasePlans_.size() != system_.phases().size()) {
        throw std::runtime_error(
            "Eulerian equation assembly requires bound AssemblyPlans.");
    }
    if (!std::isfinite(dt) || dt <= 0.0) {
        throw std::runtime_error(
            "Eulerian equation assembly requires finite positive dt.");
    }
    currentTimeStep_ = dt;
    for (size_t phase = 0; phase < system_.phases().size(); ++phase) {
        const auto& mass = system_.phases()[phase].primary.phaseMass.values();
        auto& diagonal = workspace_.momentumDiagonal[phase].values();
        for (size_t cell = 0; cell < mass.size(); ++cell) {
            diagonal[cell] = mass[cell] / dt;
            if (!std::isfinite(diagonal[cell]) || diagonal[cell] <= 0.0) {
                throw std::runtime_error(
                    "Eulerian momentum diagonal is invalid.");
            }
        }
    }
}

void PhaseEquationAssembler::assembleContinuity(
        double dt, StepSummary& summary) {
    const Field& field = system_.geometry();
    const size_t reference = system_.referencePhaseIndex();
    for (size_t phase = 0; phase < system_.phases().size(); ++phase) {
        requireTerm(*phasePlans_.at(phase).continuity,
                    Equation::TermKind::Divergence,
                    "continuity face divergence");
        if (phase == reference) continue;
        auto& mass = system_.phases()[phase].primary.phaseMass;
        const auto& old = workspace_.previousPhaseMass[phase].values();
        const auto& source = system_.sources().mass[phase];
        for (int cell : rowMap_.localCells) {
            int i = 0, j = 0, k = 0;
            field.getIJK(cell, i, j, k);
            const double rhs = -faceDivergence(
                field, workspace_.faceFlux[phase].mass, i, j, k,
                system_.config().eulerianEulerian.axisymmetric,
                system_.config().eulerianEulerian.radialCoordinate)
                + source.values()[(size_t)cell];
            const double next = old[(size_t)cell] + dt * rhs;
            if (!std::isfinite(next) || next <= 0.0) {
                throw std::runtime_error(
                    "Eulerian continuity produced non-positive phaseMass.");
            }
            // 连续方程与动量/焓矩阵分步求解。这里按质量比临时缩放主状态，
            // 仅用于连续步后的 primitive recovery 保持 U_k、h_k 不跳变；
            // 后续矩阵仍以 workspace 中的旧守恒量和完整源项装配，因此不会
            // 重复计入质量传递携带的动量或能量。
            const double oldMass = mass.values()[(size_t)cell];
            if (!std::isfinite(oldMass) || oldMass <= 0.0) {
                throw std::runtime_error(
                    "Eulerian continuity received invalid old phaseMass.");
            }
            const double ratio = next / oldMass;
            for (int component = 0; component < 3; ++component) {
                system_.phases()[phase].primary.momentum
                    [(size_t)component].values()[(size_t)cell] *= ratio;
            }
            system_.phases()[phase].primary.phaseEnthalpy
                .values()[(size_t)cell] *= ratio;
            mass.values()[(size_t)cell] = next;
            summary.massImbalance += std::abs(
                (next - old[(size_t)cell]) - dt * rhs);
        }
    }
    system_.recoverPrimitiveState();
}

void PhaseEquationAssembler::solveMomentumPredictors(double dt) {
    const Field& field = system_.geometry();
    const auto& couplings = system_.sources().momentumCouplings;
    for (size_t phase = 0; phase < system_.phases().size(); ++phase) {
        requireTerm(*phasePlans_.at(phase).momentum,
                    Equation::TermKind::Divergence,
                    "momentum convection");
        auto& state = system_.phases()[phase];
        const auto& properties = system_.phaseProperties(phase);
        ScalarField diffusivity;
        diffusivity.setupLike(field, "alphaMu." + state.name);
        for (int cell = 0; cell < field.TotalSize(); ++cell) {
            const double temperature = state.primitive.temperature
                .values()[(size_t)cell];
            diffusivity.values()[(size_t)cell] =
                state.primitive.alpha.values()[(size_t)cell]
                * Physics::Multiphase::phaseViscosity(
                    properties, temperature);
            if (turbulence_) {
                diffusivity.values()[(size_t)cell] +=
                    state.primitive.alpha.values()[(size_t)cell]
                    * turbulence_->momentumEddyViscosity(phase, cell);
            }
        }
        for (int component = 0; component < 3; ++component) {
            ScalarField source;
            source.setupLike(field, "momentumRhs." + state.name);
            for (int cell = 0; cell < field.TotalSize(); ++cell) {
                double value = system_.sources().momentum[phase]
                    [(size_t)component].values()[(size_t)cell];
                for (const auto& coupling : couplings) {
                    size_t other = system_.phases().size();
                    if (phase == coupling.first) other = coupling.second;
                    else if (phase == coupling.second) other = coupling.first;
                    if (other == system_.phases().size()) continue;
                    value -= coupling.coefficient.values()[(size_t)cell]
                        * (system_.phases()[other].primitive.velocity
                              [(size_t)component].values()[(size_t)cell]
                           - state.primitive.velocity[(size_t)component]
                                 .values()[(size_t)cell]);
                }
                source.values()[(size_t)cell] = value;
            }
            for (int cell : rowMap_.localCells) {
                int i = 0, j = 0, k = 0;
                field.getIJK(cell, i, j, k);
                const auto pressureGradient = Ops::gradient(
                    field, system_.sharedPressure(), i, j, k);
                source.values()[(size_t)cell] -=
                    state.primitive.alpha.values()[(size_t)cell]
                    * pressureGradient[(size_t)component];
            }
            auto matrix = assembleTransportEquation(
                field, rowMap_, state.primitive.velocity[(size_t)component],
                workspace_.previousMomentum[phase][(size_t)component],
                state.primary.phaseMass, diffusivity, source,
                workspace_.faceFlux[phase].mass, dt,
                system_.config().eulerianEulerian.axisymmetric,
                system_.config().eulerianEulerian.radialCoordinate,
                workspace_.momentumDiagonal[phase]);
            auto& session = momentumSolvers_[3 * phase
                                             + (size_t)component];
            LinearAlgebra::SolveResult result;
            try {
                result = session->solve(matrix);
            } catch (const std::exception& error) {
                throw std::runtime_error(
                    "Eulerian momentum solve failed for phase '"
                    + state.name + "', component="
                    + std::to_string(component) + ": " + error.what());
            }
            auto& predicted = workspace_.predictedVelocity[phase]
                                  [(size_t)component].values();
            predicted = state.primitive.velocity[(size_t)component].values();
            for (size_t row = 0; row < rowMap_.localCells.size(); ++row) {
                const int cell = rowMap_.localCells[row];
                const double old = predicted[(size_t)cell];
                predicted[(size_t)cell] = old
                    + config_.coupling.momentumRelaxation
                    * (result.solution[row] - old);
            }
        }
        for (double& value : workspace_.momentumDiagonal[phase].values()) {
            value /= config_.coupling.momentumRelaxation;
        }
    }
}

void PhaseEquationAssembler::applySemiImplicitInterphase(double dt) {
    if (!std::isfinite(dt) || dt <= 0.0) {
        throw std::runtime_error(
            "Semi-implicit interphase coupling requires positive dt.");
    }
    auto& phases = system_.phases();
    const auto& couplings = system_.sources().momentumCouplings;
    const size_t phaseCount = phases.size();
    for (int cell = 0; cell < system_.geometry().TotalSize(); ++cell) {
        std::vector<std::vector<double>> matrix(
            phaseCount, std::vector<double>(phaseCount, 0.0));
        std::vector<double> baseDiagonal(phaseCount, 0.0);
        for (size_t phase = 0; phase < phaseCount; ++phase) {
            baseDiagonal[phase] = workspace_.momentumDiagonal[phase]
                                      .values()[(size_t)cell];
            if (!std::isfinite(baseDiagonal[phase])
                || baseDiagonal[phase] <= 0.0) {
                throw std::runtime_error(
                    "Semi-implicit coupling received invalid momentum diagonal.");
            }
            matrix[phase][phase] = baseDiagonal[phase];
        }
        for (const auto& coupling : couplings) {
            const double coefficient =
                coupling.coefficient.values()[(size_t)cell];
            if (!std::isfinite(coefficient) || coefficient < 0.0) {
                throw std::runtime_error(
                    "Semi-implicit pair coefficient is invalid.");
            }
            matrix[coupling.first][coupling.first] += coefficient;
            matrix[coupling.second][coupling.second] += coefficient;
            matrix[coupling.first][coupling.second] -= coefficient;
            matrix[coupling.second][coupling.first] -= coefficient;
        }
        for (int component = 0; component < 3; ++component) {
            std::vector<double> rhs(phaseCount, 0.0);
            for (size_t phase = 0; phase < phaseCount; ++phase) {
                rhs[phase] = baseDiagonal[phase]
                    * workspace_.predictedVelocity[phase][(size_t)component]
                          .values()[(size_t)cell];
            }
            std::vector<double> solution;
            if (!Math::LinearSystem::solve(matrix, rhs, solution)) {
                throw std::runtime_error(
                    "Semi-implicit phase-pair momentum block is singular.");
            }
            for (size_t phase = 0; phase < phaseCount; ++phase) {
                if (!std::isfinite(solution[phase])) {
                    throw std::runtime_error(
                        "Semi-implicit phase-pair solution is non-finite.");
                }
                phases[phase].primitive.velocity[(size_t)component]
                    .values()[(size_t)cell] = solution[phase];
            }
        }
        for (size_t phase = 0; phase < phaseCount; ++phase) {
            workspace_.momentumDiagonal[phase].values()[(size_t)cell]
                = matrix[phase][phase];
        }
    }
    for (size_t phase = 0; phase < phases.size(); ++phase) {
        for (int cell = 0; cell < system_.geometry().TotalSize(); ++cell) {
            const double mass = phases[phase].primary.phaseMass
                                    .values()[(size_t)cell];
            for (int component = 0; component < 3; ++component) {
                const double velocity = phases[phase].primitive.velocity
                    [(size_t)component].values()[(size_t)cell];
                phases[phase].primary.momentum[(size_t)component]
                    .values()[(size_t)cell] = mass * velocity;
            }
            const double temperature = phases[phase].primitive.temperature
                                           .values()[(size_t)cell];
            const auto& properties = system_.phaseProperties(phase);
            phases[phase].primary.phaseEnthalpy.values()[(size_t)cell] =
                mass * Physics::Multiphase::phaseSpecificEnthalpy(
                    properties, temperature);
        }
    }
    system_.recoverPrimitiveState();
}

void PhaseEquationAssembler::solvePhaseEnergy(double dt) {
    const Field& field = system_.geometry();
    for (size_t phase = 0; phase < system_.phases().size(); ++phase) {
        requireTerm(*phasePlans_.at(phase).enthalpy,
                    Equation::TermKind::Divergence,
                    "enthalpy convection");
        auto& state = system_.phases()[phase];
        const auto& properties = system_.phaseProperties(phase);
        ScalarField diffusivity, source, diagonal;
        diffusivity.setupLike(field, "alphaK." + state.name);
        source.setupLike(field, "energyRhs." + state.name);
        diagonal.setupLike(field, "energyDiagonal." + state.name);
        for (int cell = 0; cell < field.TotalSize(); ++cell) {
            const double temperature = state.primitive.temperature
                .values()[(size_t)cell];
            const double cp =
                Physics::Multiphase::phaseSpecificHeat(
                    properties, temperature);
            const double conductivity =
                Physics::Multiphase::phaseThermalConductivity(
                    properties, temperature);
            diffusivity.values()[(size_t)cell] = state.primitive.alpha
                .values()[(size_t)cell]
                * conductivity / cp;
            if (turbulence_) {
                diffusivity.values()[(size_t)cell] +=
                    state.primitive.alpha.values()[(size_t)cell]
                    * turbulence_->energyEddyDiffusivity(phase, cell);
            }
            double pressureWork = 0.0;
            int i = 0, j = 0, k = 0;
            field.getIJK(cell, i, j, k);
            if (Ops::solved(
                    field, i, j, k,
                    system_.config().eulerianEulerian.axisymmetric,
                    system_.config().eulerianEulerian.radialCoordinate)) {
                const auto pressureGradient = Ops::gradient(
                    field, system_.sharedPressure(), i, j, k);
                double materialPressureRate =
                    (system_.sharedPressure().values()[(size_t)cell]
                     - workspace_.previousPressure.values()[(size_t)cell])
                    / dt;
                for (int component = 0; component < 3; ++component) {
                    materialPressureRate +=
                        state.primitive.velocity[(size_t)component]
                            .values()[(size_t)cell]
                        * pressureGradient[(size_t)component];
                }
                pressureWork = state.primitive.alpha
                    .values()[(size_t)cell] * materialPressureRate;
            }
            source.values()[(size_t)cell] =
                system_.sources().energy[phase].values()[(size_t)cell]
                + pressureWork;
        }
        auto matrix = assembleTransportEquation(
            field, rowMap_, state.primitive.enthalpy,
            workspace_.previousPhaseEnthalpy[phase],
            state.primary.phaseMass, diffusivity, source,
            workspace_.faceFlux[phase].mass, dt,
            system_.config().eulerianEulerian.axisymmetric,
            system_.config().eulerianEulerian.radialCoordinate,
            diagonal);
        LinearAlgebra::SolveResult result;
        try {
            result = energySolvers_[phase]->solve(matrix);
        } catch (const std::exception& error) {
            throw std::runtime_error(
                "Eulerian energy solve failed for phase '" + state.name
                + "': " + error.what());
        }
        for (size_t row = 0; row < rowMap_.localCells.size(); ++row) {
            state.primitive.enthalpy.values()
                [(size_t)rowMap_.localCells[row]] = result.solution[row];
        }
        for (int cell = 0; cell < field.TotalSize(); ++cell) {
            const double mass = state.primary.phaseMass.values()[(size_t)cell];
            const double enthalpy =
                state.primitive.enthalpy.values()[(size_t)cell];
            if (!std::isfinite(enthalpy) || enthalpy <= 0.0) {
                throw std::runtime_error(
                    "Eulerian enthalpy solve produced invalid state.");
            }
            state.primary.phaseEnthalpy.values()[(size_t)cell] =
                mass * enthalpy;
            state.primitive.temperature.values()[(size_t)cell] =
                Physics::Multiphase::phaseTemperatureFromEnthalpy(
                    properties, enthalpy);
        }
    }
    system_.recoverPrimitiveState();
}

} // namespace SF::EulerianEulerian
