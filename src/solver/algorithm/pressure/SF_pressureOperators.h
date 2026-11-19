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
#include "solver/algorithm/SF_faceBinding.h"

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
                      SolverAlgorithm::FaceFluxExecution& faces);

    void bind(System::StateRealization& state, FDM::IExecutionRuntime& runtime);
    double prepare(double maximumTimeStep);
    void beginFixedTimeStep();
    void beginIteration();
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
    /// @brief 可选只读验收快照，使用宿主现有 step/time；不改变状态、通量或求解顺序。
    void traceState(const char* phase,int step,double time) const;

    double iterationDelta() const { return iterationDelta_; }
    double velocityIterationDelta() const { return velocityDelta_; }
    double pressureIterationDelta() const { return pressureDelta_; }
    double fluxIterationDelta() const { return fluxDelta_; }
    double continuityDefect() const { return summary_.maxDivergenceAfter; }
    double candidateContinuityDefect() const {
        return candidateContinuityDefect_;
    }
    bool fixedTimeActive() const { return fixedTimeActive_; }

private:
    const PressureOperatorConfig config_;
    SolverAlgorithm::FaceFluxExecution& faces_;
    SolverAlgorithm::CellEquationExecution& predictor_;
    std::vector<double>& rAU_;
    const std::vector<double>& hByA_;
    FDM::IExecutionRuntime* runtime_ = nullptr;
    std::unique_ptr<LinearAlgebra::StaticDistributedNumbering> numbering_;
    std::unique_ptr<LinearAlgebra::DistributedLinearSystem> pressureSolver_;
    std::vector<std::int64_t> entityOfCell_;
    std::vector<std::int64_t> backendRow_;
    std::vector<int> candidateCells_;
    std::int64_t referenceEntity_ = -1;
    bool hasDirichlet_ = false;
    const std::array<bool,6>& physicalSides_;
    void initializeNumbering();
    void synchronize(std::vector<double>& values, int components,
                     State::ExchangeKind kind = State::ExchangeKind::State);
    void synchronize(State::DistributedFieldView& view);
    LinearAlgebra::GlobalDofId dof(int cell) const;
    void restoreGauge();
    void requireCollectively(bool valid,const char* reason) const;
    void traceAssembly() const;
    mutable unsigned long long traceAssemblyIndex_ = 0;
    mutable unsigned long long traceStateIndex_ = 0;
    State::DistributedFieldView correctionView_;
    Field* geometry_ = nullptr;
    State::DistributedFieldView* velocity_ = nullptr;
    State::DistributedFieldView* physicalVelocity_ = nullptr;
    void publishVelocity();
    State::DistributedFieldView* pressure_ = nullptr;
    double dt_ = 0.0;
    std::vector<int> cells_;
    std::vector<int> rowOfCell_;
    const std::vector<unsigned char>& boundaryCell_;
    const std::vector<unsigned char>& pressureDirichlet_;
    const std::vector<double>& iterationVelocity_;
    std::vector<double> iterationPressure_;
    std::vector<double> iterationFlux_;
    std::vector<double> correction_;
    const std::vector<double>& predictedFlux_;
    const std::vector<double>& faceResponse_;
    std::vector<double>& correctedFlux_;
    LinearAlgebra::GlobalDofSystem pressureMatrix_;
    LinearAlgebra::SolveResult pressureResult_;
    PressureOperationSummary summary_;
    bool fixedTimeActive_ = false;
    bool iterationActive_ = false;
    double iterationDelta_ = 0.0;
    double velocityDelta_ = 0.0;
    double pressureDelta_ = 0.0;
    double fluxDelta_ = 0.0;
    double candidateContinuityDefect_ = 0.0;

    void applyVelocityBoundary();
    void applyPressureBoundary();
    std::array<double,3> gradient(const std::vector<double>& values,
                                  int cell) const;
    double faceFlux(const std::vector<double>& flux, int axis,
                    int lowerCell) const;
    double divergence(const std::vector<double>& flux, int cell) const;
    double maxDivergence(const std::vector<double>& flux) const;
    bool solved(int cell) const;
    int adjacent(int cell, int axis, int sign) const;
    std::array<double,3> faceArea(int axis, int lowerCell) const;
    std::size_t faceIndex(int axis, int lowerCell) const;
};

} // namespace SF::Pressure
