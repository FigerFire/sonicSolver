/// @file SF_pressure.cpp
/// @brief 双欧拉共享压力方程及所有相压力贡献的装配。

#include "solver/equation/eulerian/SF_equations.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <stdexcept>

namespace SF::EulerianEulerian {
namespace {

struct FaceGeometry {
    int lowerI = 0, lowerJ = 0, lowerK = 0;
    int upperI = 0, upperJ = 0, upperK = 0;
    int lowerCell = -1, upperCell = -1;
    double spacing = 0.0;
    double metricNorm = 0.0;
    std::array<double, 3> area{0.0, 0.0, 0.0};
    std::array<double, 3> normal{0.0, 0.0, 0.0};
};

FaceGeometry faceGeometry(const Field& field, int axis,
                          int lowerI, int lowerJ, int lowerK,
                          bool axisymmetric, int radialCoordinate) {
    FaceGeometry face;
    face.lowerI = lowerI;
    face.lowerJ = lowerJ;
    face.lowerK = lowerK;
    int di = 0, dj = 0, dk = 0;
    Ops::offset(axis, 1, di, dj, dk);
    face.upperI = lowerI + di;
    face.upperJ = lowerJ + dj;
    face.upperK = lowerK + dk;
    face.lowerCell = field.getIdx(lowerI, lowerJ, lowerK);
    face.upperCell = field.getIdx(
        face.upperI, face.upperJ, face.upperK);
    face.spacing = Ops::spacing(
        field, lowerI, lowerJ, lowerK,
        face.upperI, face.upperJ, face.upperK);
    face.area = Ops::faceCofactor(
        field, axis, lowerI, lowerJ, lowerK,
        axisymmetric, radialCoordinate);
    for (int component = 0; component < 3; ++component) {
        const double lowerCoordinate = component == 0
            ? field.X(lowerI, lowerJ, lowerK)
            : component == 1 ? field.Y(lowerI, lowerJ, lowerK)
                             : field.Z(lowerI, lowerJ, lowerK);
        const double upperCoordinate = component == 0
            ? field.X(face.upperI, face.upperJ, face.upperK)
            : component == 1 ? field.Y(face.upperI, face.upperJ, face.upperK)
                             : field.Z(face.upperI, face.upperJ, face.upperK);
        face.normal[(size_t)component] =
            (upperCoordinate - lowerCoordinate) / face.spacing;
        face.metricNorm += face.area[(size_t)component]
                         * face.area[(size_t)component];
    }
    face.metricNorm = std::sqrt(face.metricNorm);
    if (!std::isfinite(face.metricNorm) || face.metricNorm <= 0.0) {
        throw std::runtime_error(
            "Eulerian face metric is degenerate.");
    }
    return face;
}

double pressureResponse(
        const Physics::PhaseSystems::PhaseState& phase,
        const ScalarField& diagonal, const FaceGeometry& face) {
    const double alphaLower = phase.primitive.alpha
        .values()[(size_t)face.lowerCell];
    const double alphaUpper = phase.primitive.alpha
        .values()[(size_t)face.upperCell];
    const double diagLower = diagonal.values()[(size_t)face.lowerCell];
    const double diagUpper = diagonal.values()[(size_t)face.upperCell];
    if (!std::isfinite(diagLower) || diagLower <= 0.0
        || !std::isfinite(diagUpper) || diagUpper <= 0.0) {
        throw std::runtime_error(
            "Rhie-Chow interpolation found invalid momentum diagonal.");
    }
    const double response = 0.5
        * (alphaLower * alphaLower / diagLower
           + alphaUpper * alphaUpper / diagUpper);
    if (!std::isfinite(response) || response <= 0.0) {
        throw std::runtime_error(
            "Rhie-Chow pressure response is invalid.");
    }
    return response;
}

void visitAdjacentFaces(const Field& field,
                        const std::function<void(int,int,int,int)>& visitor,
                        bool axisymmetric, int radialCoordinate) {
    const int ng = field.NG();
    for (int k = ng; k < ng + field.NZ(); ++k) {
        for (int j = ng; j < ng + field.NY(); ++j) {
            for (int i = ng; i < ng + field.NX(); ++i) {
                if (!Ops::solved(
                        field, i, j, k,
                        axisymmetric, radialCoordinate)) continue;
                for (int axis = 0; axis < 3; ++axis) {
                    if (!Math::isDirectionActiveIndex(axis)) continue;
                    visitor(axis, i, j, k);
                    int di = 0, dj = 0, dk = 0;
                    Ops::offset(axis, -1, di, dj, dk);
                    visitor(axis, i + di, j + dj, k + dk);
                }
            }
        }
    }
}

double divergence(const Field& field, const std::vector<double>& flux,
                  int i, int j, int k,
                  bool axisymmetric, int radialCoordinate) {
    double result = 0.0;
    for (int axis = 0; axis < 3; ++axis) {
        if (!Math::isDirectionActiveIndex(axis)) continue;
        int di = 0, dj = 0, dk = 0;
        Ops::offset(axis, -1, di, dj, dk);
        result += flux[PhaseFaceFlux::index(field, axis, i, j, k)]
                - flux[PhaseFaceFlux::index(
                    field, axis, i + di, j + dj, k + dk)];
    }
    return Ops::inverseCellVolume(
        field, i, j, k, axisymmetric, radialCoordinate) * result;
}

} // namespace

void PhaseEquationAssembler::buildMomentumInterpolatedFlux() {
    const Field& field = system_.geometry();
    const ScalarField& pressure = system_.sharedPressure();
    const auto& geometry = system_.config().eulerianEulerian;
    for (size_t phaseIndex = 0;
         phaseIndex < system_.phases().size(); ++phaseIndex) {
        const auto& phase = system_.phases()[phaseIndex];
        auto& flux = workspace_.faceFlux[phaseIndex];
        std::fill(flux.volume.begin(), flux.volume.end(), 0.0);
        std::fill(flux.mass.begin(), flux.mass.end(), 0.0);
        visitAdjacentFaces(field, [&](int axis, int i, int j, int k) {
            const FaceGeometry face = faceGeometry(
                field, axis, i, j, k,
                geometry.axisymmetric, geometry.radialCoordinate);
            const size_t index = PhaseFaceFlux::index(
                field, axis, i, j, k);
            const double alpha = 0.5
                * (phase.primitive.alpha.values()[(size_t)face.lowerCell]
                   + phase.primitive.alpha.values()[(size_t)face.upperCell]);
            const double density = 0.5
                * (phase.primitive.density.values()[(size_t)face.lowerCell]
                   + phase.primitive.density.values()[(size_t)face.upperCell]);
            double velocityFlux = 0.0;
            for (int component = 0; component < 3; ++component) {
                const double velocity = 0.5
                    * (phase.primitive.velocity[(size_t)component]
                           .values()[(size_t)face.lowerCell]
                       + phase.primitive.velocity[(size_t)component]
                           .values()[(size_t)face.upperCell]);
                velocityFlux += velocity * face.area[(size_t)component];
            }
            const double directGradient =
                (pressure.values()[(size_t)face.upperCell]
                 - pressure.values()[(size_t)face.lowerCell]) / face.spacing;
            double interpolatedGradient = directGradient;
            if (Ops::solved(
                    field, face.lowerI, face.lowerJ, face.lowerK,
                    geometry.axisymmetric, geometry.radialCoordinate)
                && Ops::solved(
                    field, face.upperI, face.upperJ, face.upperK,
                    geometry.axisymmetric, geometry.radialCoordinate)) {
                const auto lowerGradient = Ops::gradient(
                    field, pressure, face.lowerI, face.lowerJ, face.lowerK);
                const auto upperGradient = Ops::gradient(
                    field, pressure, face.upperI, face.upperJ, face.upperK);
                interpolatedGradient = 0.0;
                for (int component = 0; component < 3; ++component) {
                    interpolatedGradient += 0.5
                        * (lowerGradient[(size_t)component]
                           + upperGradient[(size_t)component])
                        * face.normal[(size_t)component];
                }
            }
            const double response = pressureResponse(
                phase, workspace_.momentumDiagonal[phaseIndex], face);
            const double volumeFlux = alpha * velocityFlux
                - response * (directGradient - interpolatedGradient)
                    * face.metricNorm;
            if (!std::isfinite(volumeFlux)
                || !std::isfinite(density) || density <= 0.0) {
                throw std::runtime_error(
                    "Rhie-Chow interpolation produced invalid phase flux.");
            }
            flux.volume[index] = volumeFlux;
            flux.mass[index] = density * volumeFlux;
        }, geometry.axisymmetric, geometry.radialCoordinate);
    }
}

LinearAlgebra::SolveResult
PhaseEquationAssembler::solvePressureCorrection() {
    if (!pressurePlan_) {
        throw std::runtime_error(
            "Eulerian pressure correction requires a bound AssemblyPlan.");
    }
    requireTerm(*pressurePlan_,Equation::TermKind::Constraint,
                "shared-pressure correction");
    const Field& field = system_.geometry();
    const auto& geometry = system_.config().eulerianEulerian;
    if (rowMap_.total <= 0) {
        throw std::runtime_error(
            "Eulerian pressure equation has no unknowns.");
    }
    if (config_.referenceCell < 0
        || config_.referenceCell >= rowMap_.total) {
        throw std::runtime_error(
            "solverProperties referenceCell is outside pressure matrix.");
    }
    LinearAlgebra::SparseSystem matrix;
    matrix.firstRow = rowMap_.first;
    matrix.lastRow = rowMap_.last;
    matrix.globalSize = rowMap_.total;
    matrix.rows.reserve(rowMap_.localCells.size());
    matrix.rhs.reserve(rowMap_.localCells.size());
    matrix.initialGuess.assign(rowMap_.localCells.size(), 0.0);
    for (int cell : rowMap_.localCells) {
        int i = 0, j = 0, k = 0;
        field.getIJK(cell, i, j, k);
        LinearAlgebra::SparseRow row;
        row.globalRow = rowMap_.global[(size_t)cell];
        if (row.globalRow == config_.referenceCell) {
            row.columns = {row.globalRow};
            row.values = {1.0};
            matrix.rows.push_back(std::move(row));
            matrix.rhs.push_back(0.0);
            continue;
        }
        double diagonal = 0.0;
        for (int axis = 0; axis < 3; ++axis) {
            if (!Math::isDirectionActiveIndex(axis)) continue;
            for (int sign : {-1, 1}) {
                int di = 0, dj = 0, dk = 0;
                Ops::offset(axis, sign, di, dj, dk);
                const int ni = i + di, nj = j + dj, nk = k + dk;
                const int lowerI = sign > 0 ? i : ni;
                const int lowerJ = sign > 0 ? j : nj;
                const int lowerK = sign > 0 ? k : nk;
                const FaceGeometry face = faceGeometry(
                    field, axis, lowerI, lowerJ, lowerK,
                    geometry.axisymmetric, geometry.radialCoordinate);
                double response = 0.0;
                for (size_t phase = 0;
                     phase < system_.phases().size(); ++phase) {
                    response += pressureResponse(
                        system_.phases()[phase],
                        workspace_.momentumDiagonal[phase], face);
                }
                const double coefficient =
                    Ops::inverseCellVolume(
                        field, i, j, k,
                        geometry.axisymmetric, geometry.radialCoordinate)
                    * response * face.metricNorm / face.spacing;
                if (!std::isfinite(coefficient) || coefficient <= 0.0) {
                    throw std::runtime_error(
                        "Eulerian pressure coefficient is invalid.");
                }
                diagonal += coefficient;
                const int neighbour = field.getIdx(ni, nj, nk);
                const auto global = rowMap_.global[(size_t)neighbour];
                if (global >= 0) {
                    row.columns.push_back(global);
                    row.values.push_back(-coefficient);
                }
            }
        }
        row.columns.push_back(row.globalRow);
        row.values.push_back(diagonal);
        double mixtureDivergence = 0.0;
        for (const auto& flux : workspace_.faceFlux) {
            mixtureDivergence += divergence(
                field, flux.volume, i, j, k,
                geometry.axisymmetric, geometry.radialCoordinate);
        }
        double volumeSource = 0.0;
        double compressibility = 0.0;
        for (size_t phase = 0;
             phase < system_.phases().size(); ++phase) {
            const auto& state = system_.phases()[phase];
            const double density =
                state.primitive.density.values()[(size_t)cell];
            const double alpha =
                state.primitive.alpha.values()[(size_t)cell];
            volumeSource += system_.sources().mass[phase]
                .values()[(size_t)cell] / density;
            const auto& properties = system_.phaseProperties(phase);
            if (Physics::Multiphase::normalizeModelType(
                    properties.thermoModel) == "perfectgas"
                && Physics::Multiphase::normalizeModelType(
                    properties.densityLaw.model)
                    != "polynomialtemperature") {
                const double pressure = system_.sharedPressure()
                    .values()[(size_t)cell];
                compressibility += alpha / pressure;
            }
        }
        const double pressureHistory = compressibility
            * (system_.sharedPressure().values()[(size_t)cell]
               - workspace_.previousPressure.values()[(size_t)cell])
            / currentTimeStep_;
        diagonal += compressibility / currentTimeStep_;
        row.values.back() = diagonal;
        matrix.rows.push_back(std::move(row));
        matrix.rhs.push_back(
            volumeSource - mixtureDivergence - pressureHistory);
    }
    return pressureSolver_.solve(matrix);
}

void PhaseEquationAssembler::setPressureCorrection(
        const std::vector<double>& correction) {
    if (correction.size() != rowMap_.localCells.size()) {
        throw std::runtime_error(
            "Pressure correction does not match local matrix rows.");
    }
    std::fill(workspace_.pressureCorrection.values().begin(),
              workspace_.pressureCorrection.values().end(), 0.0);
    for (size_t row = 0; row < correction.size(); ++row) {
        workspace_.pressureCorrection.values()
            [(size_t)rowMap_.localCells[row]] = correction[row];
    }
}

void PhaseEquationAssembler::correctAllPhases() {
    const Field& field = system_.geometry();
    const ScalarField& correction = workspace_.pressureCorrection;
    for (int cell : rowMap_.localCells) {
        int i = 0, j = 0, k = 0;
        field.getIJK(cell, i, j, k);
        const auto gradient = Ops::gradient(field, correction, i, j, k);
        for (size_t phase = 0; phase < system_.phases().size(); ++phase) {
            auto& state = system_.phases()[phase];
            const double alpha = state.primitive.alpha
                                     .values()[(size_t)cell];
            const double diagonal = workspace_.momentumDiagonal[phase]
                                        .values()[(size_t)cell];
            const double mass = state.primary.phaseMass
                                    .values()[(size_t)cell];
            for (int component = 0; component < 3; ++component) {
                const double velocity = state.primitive.velocity
                    [(size_t)component].values()[(size_t)cell]
                    - alpha / diagonal * gradient[(size_t)component];
                if (!std::isfinite(velocity)) {
                    throw std::runtime_error(
                        "Pressure correction produced invalid phase velocity.");
                }
                state.primary.momentum[(size_t)component]
                    .values()[(size_t)cell] = mass * velocity;
            }
        }
        const double pressure = system_.sharedPressure()
                                    .values()[(size_t)cell]
            + config_.pressureRelaxation
                * correction.values()[(size_t)cell];
        if (!std::isfinite(pressure) || pressure <= 0.0) {
            throw std::runtime_error(
                "Pressure correction produced non-positive pressure.");
        }
        system_.sharedPressure().values()[(size_t)cell] = pressure;
    }
    system_.recoverPrimitiveState();
}

void PhaseEquationAssembler::correctCanonicalFaceFlux() {
    const Field& field = system_.geometry();
    const ScalarField& correction = workspace_.pressureCorrection;
    const auto& geometry = system_.config().eulerianEulerian;
    for (size_t phase = 0; phase < system_.phases().size(); ++phase) {
        visitAdjacentFaces(field, [&](int axis, int i, int j, int k) {
            const FaceGeometry face = faceGeometry(
                field, axis, i, j, k,
                geometry.axisymmetric, geometry.radialCoordinate);
            const size_t index = PhaseFaceFlux::index(
                field, axis, i, j, k);
            const double response = pressureResponse(
                system_.phases()[phase],
                workspace_.momentumDiagonal[phase], face);
            const double pressureGradient =
                (correction.values()[(size_t)face.upperCell]
                 - correction.values()[(size_t)face.lowerCell])
                / face.spacing;
            const double deltaFlux =
                response * pressureGradient * face.metricNorm;
            workspace_.faceFlux[phase].volume[index] -= deltaFlux;
            const double density = 0.5
                * (system_.phases()[phase].primitive.density
                       .values()[(size_t)face.lowerCell]
                   + system_.phases()[phase].primitive.density
                       .values()[(size_t)face.upperCell]);
            workspace_.faceFlux[phase].mass[index] -= density * deltaFlux;
        }, geometry.axisymmetric, geometry.radialCoordinate);
    }
}

} // namespace SF::EulerianEulerian
