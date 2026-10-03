#pragma once

/// @file SF_pressureOperators.h
/// @brief Constant-density pressure-constraint numerical leaves.
///
/// Each method implements one CompiledSolvePlan operation. This class owns no
/// corrector or timestep loop; PlanExecutor alone repeats the operations.

#include "solver/system/SF_numericalSystem.h"
#include "core/system/SF_solveProgram.h"
#include "solver/system/SF_stateRealizer.h"
#include "solver/linearAlgebra/SF_globalDofSystem.h"
#include "core/config/types/SF_pressureConfigTypes.h"
#include "core/state/SF_valueTypes.h"

#include "core/interfaces/SF_executionRuntime.h"

#include <array>
#include <memory>
#include <vector>

namespace SF::Pressure {

using PressureOperatorConfig = System::PressureNumericalConfig;

struct PressureOperationSummary {
    double maxDivergenceBefore = 0.0;
    double maxDivergenceAfter = 0.0;
    int iterations = 0;
    double relativeResidual = 0.0;
};

class PressureOperators {
public:
    PressureOperators(PressureOperatorConfig config,
                      const System::CompiledNumericalSystem& numerics);

    void bind(System::StateRealization& state, FDM::IExecutionRuntime& runtime);
    double prepare(double maximumTimeStep);
    void beginFixedTimeStep();
    void beginIteration();
    void assembleMomentum();
    void solveMomentum(System::TargetKind target = System::TargetKind::Physical);
    void preparePressureBoundary();
    void assemblePressure();
    void solvePressure();
    void preparePressureUpdate();
    void correctVelocity();
    void correctFlux();
    PressureOperationSummary commitCorrection();
    void applyRelaxation();
    void restoreFluxConsistency();
    void evaluateConvergence();
    void endIteration();
    void commitStep();

    const std::vector<double>& correctedFaceFlux() const {
        return correctedFlux_;
    }
    double iterationDelta() const { return iterationDelta_; }
    double velocityIterationDelta() const { return velocityDelta_; }
    double pressureIterationDelta() const { return pressureDelta_; }
    double fluxIterationDelta() const { return fluxDelta_; }
    double continuityDefect() const { return summary_.maxDivergenceAfter; }
    double candidateContinuityDefect() const {
        return candidateContinuityDefect_;
    }
    bool converged() const { return converged_; }
    bool fixedTimeActive() const { return fixedTimeActive_; }

private:
    const PressureOperatorConfig config_;
    const System::CompiledNumericalSystem& numerics_;
    FDM::IExecutionRuntime* runtime_ = nullptr;
    std::unique_ptr<LinearAlgebra::StaticDistributedNumbering> numbering_;
    std::unique_ptr<LinearAlgebra::DistributedLinearSystem> pressureSolver_;
    std::vector<std::int64_t> entityOfCell_;
    std::vector<std::int64_t> backendRow_;
    std::vector<System::CompiledSpatialBinding> momentumTerms_;
    std::vector<int> candidateCells_;
    std::int64_t referenceEntity_ = -1;
    bool hasDirichlet_ = false;
    std::array<bool,6> physicalSides_{};
    void initializeNumbering();
    void synchronize(std::vector<double>& values, int components,
                     State::ExchangeKind kind = State::ExchangeKind::State);
    void synchronize(State::DistributedFieldView& view);
    LinearAlgebra::GlobalDofId dof(int cell) const;
    void restoreGauge();
    void requireCollectively(bool valid,const char* reason) const;
    void traceAssembly() const;
    mutable unsigned long long traceAssemblyIndex_ = 0;
    System::StateRealization* stateViews_ = nullptr;
    State::DistributedFieldView correctionView_;
    State::DistributedFieldView fluxView_;
    Field* geometry_ = nullptr;
    State::DistributedFieldView* velocity_ = nullptr;
    State::DistributedFieldView* physicalVelocity_ = nullptr;
    State::DistributedFieldView workingVelocityView_;
    std::vector<double> workingVelocity_;
    void publishVelocity();
    State::DistributedFieldView* pressure_ = nullptr;
    double density_ = 0.0;
    double dt_ = 0.0;
    std::vector<int> cells_;
    std::vector<int> rowOfCell_;
    std::vector<unsigned char> boundaryCell_;
    std::vector<unsigned char> pressureDirichlet_;
    std::vector<double> hByA_;
    std::vector<double> rAU_;
    std::vector<double> baseVelocity_;
    std::vector<double> iterationVelocity_;
    std::vector<double> iterationPressure_;
    std::vector<double> iterationFlux_;
    std::vector<double> correction_;
    std::vector<double> predictedFlux_;
    std::vector<double> faceResponse_;
    std::vector<double> correctedFlux_;
    LinearAlgebra::GlobalDofSystem pressureMatrix_;
    LinearAlgebra::SolveResult pressureResult_;
    PressureOperationSummary summary_;
    bool fixedTimeActive_ = false;
    bool iterationActive_ = false;
    int iterationIndex_ = 0;
    double iterationDelta_ = 0.0;
    double velocityDelta_ = 0.0;
    double pressureDelta_ = 0.0;
    double fluxDelta_ = 0.0;
    double candidateContinuityDefect_ = 0.0;
    bool converged_ = false;

    void applyVelocityBoundary();
    void applyPressureBoundary();
    void reconstructFaceFlux(std::vector<double>& flux);
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

} // namespace SF::Pressure
