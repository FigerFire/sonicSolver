#pragma once

/// @file SF_pressureOperators.h
/// @brief Constant-density pressure-constraint numerical leaves.
///
/// Each method implements one CompiledSolvePlan operation. This class owns no
/// corrector or timestep loop; PlanExecutor alone repeats the operations.

#include "solver/system/SF_numericalSystem.h"
#include "solver/system/SF_stateRealizer.h"
#include "solver/linearAlgebra/SF_linearAlgebra.h"
#include "core/config/SF_configTypes.h"

#include <array>
#include <vector>

namespace SF::PressureBased {

struct PressureOperationSummary {
    double maxDivergenceBefore = 0.0;
    double maxDivergenceAfter = 0.0;
    int iterations = 0;
    double relativeResidual = 0.0;
};

class PressureOperators {
public:
    PressureOperators(const FDM::SolverConfig& config,
                      const System::CompiledNumericalSystem& numerics);

    void bind(const System::StateRealization& state);
    double prepare(double maximumTimeStep, double currentTime);
    void assembleMomentum();
    void solveMomentum();
    void preparePressureBoundary();
    void assemblePressure();
    void solvePressure();
    void preparePressureUpdate();
    void correctVelocity();
    void correctFlux();
    PressureOperationSummary commitCorrection();
    void commitStep();

    const std::vector<double>& correctedFaceFlux() const {
        return correctedFlux_;
    }

private:
    const FDM::SolverConfig& config_;
    const System::CompiledNumericalSystem& numerics_;
    LinearAlgebra::SolverSession pressureSolver_;
    Field* geometry_ = nullptr;
    State::DistributedFieldView* velocity_ = nullptr;
    State::DistributedFieldView* pressure_ = nullptr;
    double density_ = 0.0;
    double dt_ = 0.0;
    std::vector<int> cells_;
    std::vector<int> rowOfCell_;
    std::vector<unsigned char> boundaryCell_;
    std::vector<unsigned char> pressureDirichlet_;
    std::vector<double> hByA_;
    std::vector<double> rAU_;
    std::vector<double> correction_;
    std::vector<double> predictedFlux_;
    std::vector<double> faceResponse_;
    std::vector<double> correctedFlux_;
    LinearAlgebra::SparseSystem pressureMatrix_;
    LinearAlgebra::SolveResult pressureResult_;
    PressureOperationSummary summary_;

    void applyVelocityBoundary();
    void applyPressureBoundary();
    void extendHalos(State::DistributedFieldView& view, int offset,
                     int components);
    std::array<double,3> gradient(const std::vector<double>& values,
                                  int cell) const;
    std::array<double,3> pressureGradient(int cell) const;
    double faceFlux(const std::vector<double>& flux, int axis,
                    int lowerCell) const;
    double divergence(const std::vector<double>& flux, int cell) const;
    double maxDivergence(const std::vector<double>& flux) const;
    bool solved(int cell) const;
    int adjacent(int cell, int axis, int sign) const;
    std::array<double,3> faceArea(int axis, int lowerCell) const;
    double spacing(int axis, int lowerCell) const;
    std::size_t faceIndex(int axis, int lowerCell) const;
};

} // namespace SF::PressureBased
