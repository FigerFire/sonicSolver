/// @file SF_transport.cpp
/// @brief 双欧拉面输运系数、梯度和通量的结构网格算子。

#include "solver/algorithm/eulerian/equation/SF_transport.h"
#include "SF_canonicalFace.h"
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include "solver/boundary/SF_scalarBoundary.h"
#include "solver/algorithm/eulerian/equation/SF_equations.h"

#include <cmath>
#include <stdexcept>

namespace SF::EulerianEulerian {

LinearAlgebra::SparseSystem assembleTransportEquation(
        const Field& field,
        const LinearAlgebra::DistributedRowMap& rowMap,
        const ScalarField& unknown,
        const ScalarField& previousConserved,
        const ScalarField& capacity,
        const ScalarField& diffusivity,
        const ScalarField& explicitSource,
        const std::vector<double>& canonicalMassFlux,
        double timeStep,
        bool axisymmetric,
        int radialCoordinate,
        ScalarField& diagonal,
        const ScalarField* implicitSink) {
    if(!std::isfinite(timeStep) || timeStep<=0)
        throw std::runtime_error("Eulerian transport requires finite positive dt.");
    const auto size=static_cast<size_t>(field.TotalSize());
    for(const ScalarField* value:{&unknown,&previousConserved,&capacity,&diffusivity,&explicitSource,static_cast<const ScalarField*>(&diagonal)})
        if(value->values().size()!=size)
            throw std::runtime_error("Missing or mismatched Eulerian transport storage: "+value->name());
    if((implicitSink && implicitSink->values().size()!=size) || canonicalMassFlux.size()!=3*size
        || unknown.values().data()==previousConserved.values().data())
        throw std::runtime_error("Eulerian transport old-time/face/sink layout or writable alias is invalid.");
    LinearAlgebra::SparseSystem matrix;
    matrix.firstRow = rowMap.first;
    matrix.lastRow = rowMap.last;
    matrix.globalSize = rowMap.total;
    matrix.rows.reserve(rowMap.localCells.size());
    matrix.rhs.reserve(rowMap.localCells.size());
    matrix.initialGuess.reserve(rowMap.localCells.size());

    for (int cell : rowMap.localCells) {
        int i = 0, j = 0, k = 0;
        field.getIJK(cell, i, j, k);
        const double cellCapacity = capacity.values()[(size_t)cell];
        if (!std::isfinite(cellCapacity) || cellCapacity <= 0.0) {
            throw std::runtime_error(
                "Eulerian transport capacity must remain positive.");
        }

        LinearAlgebra::SparseRow row;
        row.globalRow = rowMap.global[(size_t)cell];
        const double inverseVolume = Ops::inverseCellVolume(
            field, i, j, k, axisymmetric, radialCoordinate);
        double diagonalValue = cellCapacity / timeStep;
        if (implicitSink) {
            const double sink = implicitSink->values()[(size_t)cell];
            if (!std::isfinite(sink) || sink < 0.0) {
                throw std::runtime_error(
                    "Eulerian implicit transport sink is invalid.");
            }
            diagonalValue += sink;
        }
        double rhs = previousConserved.values()[(size_t)cell] / timeStep
                   + explicitSource.values()[(size_t)cell];

        auto addColumn = [&](std::int64_t column, double value) {
            for (size_t entry = 0; entry < row.columns.size(); ++entry) {
                if (row.columns[entry] == column) {
                    row.values[entry] += value;
                    return;
                }
            }
            row.columns.push_back(column);
            row.values.push_back(value);
        };

        for (int axis = 0; axis < 3; ++axis) {
            if (!Math::isDirectionActiveIndex(axis)) continue;
            for (int sign : {-1, 1}) {
                int di = 0, dj = 0, dk = 0;
                Ops::offset(axis, sign, di, dj, dk);
                const int ni = i + di;
                const int nj = j + dj;
                const int nk = k + dk;
                const int neighbour = field.getIdx(ni, nj, nk);
                const double spacing = Ops::spacing(
                    field, i, j, k, ni, nj, nk);
                const int lowerI = sign > 0 ? i : ni;
                const int lowerJ = sign > 0 ? j : nj;
                const int lowerK = sign > 0 ? k : nk;
                const auto area = Ops::faceCofactor(
                    field, axis, lowerI, lowerJ, lowerK,
                    axisymmetric, radialCoordinate);
                const double areaNorm = std::sqrt(
                    area[0]*area[0] + area[1]*area[1]
                    + area[2]*area[2]);
                const double gamma = 0.5
                    * (diffusivity.values()[(size_t)cell]
                       + diffusivity.values()[(size_t)neighbour]);
                const double diffusion =
                    inverseVolume * gamma * areaNorm / spacing;
                if (!std::isfinite(diffusion) || diffusion < 0.0) {
                    throw std::runtime_error(
                        "Eulerian diffusion coefficient is invalid.");
                }
                diagonalValue += diffusion;
                const auto global = rowMap.global[(size_t)neighbour];
                if (global >= 0) {
                    addColumn(global, -diffusion);
                } else {
                    rhs += diffusion
                        * unknown.values()[(size_t)neighbour];
                }

                const size_t face = PhaseFaceFlux::index(
                    field, axis, lowerI, lowerJ, lowerK);
                const double outwardMassFlux =
                    Numerics::CanonicalFace::outward(
                        canonicalMassFlux[face],
                        sign>0
                            ?Numerics::CanonicalFace::CellSide::LowerOwner
                            :Numerics::CanonicalFace::CellSide::UpperNeighbour);
                const double convection =
                    inverseVolume * outwardMassFlux;
                if (!std::isfinite(convection)) {
                    throw std::runtime_error(
                        "Eulerian canonical mass flux is non-finite.");
                }
                if (convection >= 0.0) {
                    diagonalValue += convection;
                } else if (global >= 0) {
                    addColumn(global, convection);
                } else {
                    rhs -= convection
                        * unknown.values()[(size_t)neighbour];
                }
            }
        }

        diagonal.values()[(size_t)cell] = diagonalValue;
        addColumn(row.globalRow, diagonalValue);
        matrix.rows.push_back(std::move(row));
        matrix.rhs.push_back(rhs);
        matrix.initialGuess.push_back(
            unknown.values()[(size_t)cell]);
    }
    return matrix;
}

} // namespace SF::EulerianEulerian




namespace SF::EulerianEulerian {
namespace {
std::ofstream& turbulenceEvidence() {
    static std::ofstream output;
    if (!output.is_open()) if (const char* path=std::getenv("SF_EULERIAN_TURBULENCE_EVIDENCE")) {
        output.open(path); output << std::setprecision(17);
        if (!output) throw std::runtime_error("Cannot open Eulerian turbulence evidence file.");
    }
    return output;
}
void recordField(const char* role,const ScalarField& field) {
    auto& out=turbulenceEvidence(); if (!out.is_open()) return;
    out << role << ' ' << field.name();
    for (double value:field.values()) out << ' ' << value;
    out << '\n';
}
void recordMatrix(const LinearAlgebra::DistributedRowMap& map,
        const LinearAlgebra::SparseSystem& matrix,const ScalarField& old,
        const ScalarField& diagonal) {
    auto& out=turbulenceEvidence(); if (!out.is_open()) return;
    recordField("old",old);recordField("diagonal",diagonal);
    for (size_t row=0;row<matrix.rows.size();++row) {
        const auto& value=matrix.rows[row];
        out << "row " << map.localCells[row] << ' ' << value.globalRow;
        for (size_t column=0;column<value.columns.size();++column)
            out << ' ' << value.columns[column] << ':' << value.values[column];
        out << " rhs " << matrix.rhs[row] << " initial " << matrix.initialGuess[row] << '\n';
    }
}
void recordResult(const LinearAlgebra::SolveResult& result) {
    auto& out=turbulenceEvidence();if(!out.is_open()) return;
    out << "solve " << result.iterations << ' ' << result.relativeResidual;
    for(double value:result.solution) out << ' ' << value;
    out << '\n';
}
}


void PhaseEquationAssembler::attachTurbulence(
        Turbulence::EquationSystem* turbulence) {
    turbulence_ = turbulence;
    turbulenceSolvers_.clear();
    workspace_.turbulenceOldTime.clear();
    if (!turbulence_) return;
    applyTurbulenceBoundary();
    turbulence_->prepare(system_);
    turbulence_->validate(system_,"initialization");
    if (!turbulence_->hasTransportEquations()) return;
    workspace_.turbulenceDiagonal.setupLike(system_.geometry(),"turbulenceDiagonal");
    for (auto& state : turbulence_->states()) {
        for (const auto* equation :
             turbulence_->transportEquations(state)) {
            TurbulenceTimeLevel old{state.phaseIndex,&equation->variable,{}};
            old.previousConserved.setupLike(system_.geometry(),"previous."+equation->variable.name());
            const auto& mass=system_.phases().at(state.phaseIndex).primary.phaseMass.values();
            for(size_t cell=0;cell<mass.size();++cell)
                old.previousConserved.values()[cell]=mass[cell]*equation->variable.values()[cell];
            workspace_.turbulenceOldTime.push_back(std::move(old));
            turbulenceSolvers_.push_back(
                std::make_unique<LinearAlgebra::SolverSession>(
                    config_.linear.turbulence));
        }
    }
}

void PhaseEquationAssembler::applyTurbulenceBoundary() {
    if (!turbulence_ || !turbulence_->hasTransportEquations()) return;
    const auto& config=turbulence_->config();
    for (auto& state : turbulence_->states()) {
        Boundary::applyNodeScalarBoundary(
            system_.geometry(), config.scalars.kBoundary,
            state.kineticEnergy.variable);
        ScalarField* second = nullptr;
        const std::vector<BCSetting<double>>* secondBoundary = nullptr;
        if (config.model == FDM::TurbulenceModelKind::kEpsilon) {
            second = &state.dissipation.variable;
            secondBoundary = &config.scalars.epsilonBoundary;
        } else {
            second = &state.specificDissipation.variable;
            secondBoundary = &config.scalars.omegaBoundary;
        }
        Boundary::applyNodeScalarBoundary(
            system_.geometry(), *secondBoundary, *second);
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

void PhaseEquationAssembler::beginTurbulenceStep(std::uint64_t step) {
    if(turbulenceStep_) throw std::runtime_error("Eulerian turbulence has an uncommitted physical step.");
    turbulenceStep_=step;coefficientEpoch_.reset();solvedTurbulenceOuter_.reset();
}
void PhaseEquationAssembler::publishTurbulenceCoefficients(std::uint64_t step,int outer) {
    if(turbulenceStep_!=step || outer<0)
        throw std::runtime_error("Wrong Eulerian turbulence coefficient physical epoch.");
    coefficientEpoch_={{step,outer}};
}
void PhaseEquationAssembler::commitTimeLevel(std::uint64_t step) {
    if(turbulenceStep_!=step || (turbulence_ && turbulence_->hasTransportEquations() && !solvedTurbulenceOuter_))
        throw std::runtime_error("Missing or duplicate Eulerian physical commit.");
    workspace_.commitTimeLevel(system_);
    turbulenceStep_.reset();coefficientEpoch_.reset();solvedTurbulenceOuter_.reset();
}
void PhaseEquationAssembler::solveTurbulence(double timeStep,std::uint64_t step,int outer) {
    if (!turbulence_ || !turbulence_->hasTransportEquations()) return;
    if(turbulenceStep_!=step || coefficientEpoch_!=std::make_pair(step,outer)
        || solvedTurbulenceOuter_==outer)
        throw std::runtime_error("Eulerian turbulence requires freshly published frozen coefficients for this physical step/outer; duplicate transport is invalid.");
    const Field& field = system_.geometry();
    auto& diagonal=workspace_.turbulenceDiagonal;
    size_t solverIndex = 0;

    for (auto& phaseState : turbulence_->states()) {
        const size_t phase = phaseState.phaseIndex;
        auto solveEquation = [&](
                Turbulence::TransportEquationState& equation) {
            const auto& old=workspace_.turbulenceOldTime.at(solverIndex);
            if (old.phaseIndex!=phase || old.physical!=&equation.variable || old.previousConserved.empty()
                || old.physical->values().data()==old.previousConserved.values().data())
                throw std::runtime_error("Eulerian turbulence old-time target/phase mismatch.");
            auto matrix = assembleTransportEquation(
                field, rowMap_, equation.variable,
                old.previousConserved,
                system_.phases()[phase].primary.phaseMass,
                equation.diffusivity, equation.explicitSource,
                workspace_.faceFlux[phase].mass, timeStep,
                system_.config().eulerianEulerian.axisymmetric,
                system_.config().eulerianEulerian.radialCoordinate,
                diagonal, &equation.implicitSink);
            recordField("target",equation.variable);
            recordField("diffusivity",equation.diffusivity);
            recordField("source",equation.explicitSource);
            recordField("sink",equation.implicitSink);
            recordMatrix(rowMap_,matrix,old.previousConserved,diagonal);
            LinearAlgebra::SolveResult result;
            try {
                result =
                    turbulenceSolvers_.at(solverIndex)->solve(matrix);
            } catch (const std::exception& error) {
                throw std::runtime_error(
                    "Eulerian turbulence solve failed for '"
                    + equation.variable.name() + "': " + error.what());
            }
            recordResult(result);
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

    solvedTurbulenceOuter_=outer;
}

void PhaseEquationAssembler::refreshTurbulenceClosure() {
    turbulence_->prepare(system_);
    for(auto& state:turbulence_->states()) {
        for(auto* equation:turbulence_->transportEquations(state)) recordField("published",equation->variable);
        recordField("mut",state.eddyViscosity);
    }
}

} // namespace SF::EulerianEulerian
