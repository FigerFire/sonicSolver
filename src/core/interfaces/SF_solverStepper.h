#pragma once

/// @file SF_solverStepper.h
/// @brief Stable timestep state, services, result, and stepper port.

#include "SF_boundaryPipeline.h"
#include "SF_equationCoupling.h"
#include "SF_executionRuntime.h"
#include "SF_observer.h"
#include "SF_transportModel.h"
#include "core/state/SF_stateBundle.h"

#include <limits>
#include <string>

namespace SF {
namespace System { struct CompiledSolvePlan; }
namespace FDM {

struct SolverState {
    State::StateBundle* bundle = nullptr;
    double maximumTimeStep = std::numeric_limits<double>::max();
};

struct SolverServices {
    IBoundaryPipeline* boundaryPipeline = nullptr;
    ImmersedCouplingPorts immersed;
    IExecutionRuntime* executionRuntime = nullptr;
    ITransportModel* transportModel = nullptr;
    IEquationSystemCoupling* equationSystem = nullptr;
    ISolverObserver* observer = nullptr;
};

struct StepResult {
    bool accepted = false;
    double dt = 0.0;
    double time = 0.0;
    int step = 0;
    bool outputDue = false;
    bool finished = false;
    std::string message;
};

class INavierStokesStepper {
public:
    virtual ~INavierStokesStepper() = default;
    virtual void bindServices(SolverServices services) = 0;
    virtual void bindSolvePlan(const System::CompiledSolvePlan& plan) = 0;
    virtual void prepare(SolverState& state) = 0;
    virtual StepResult advance(SolverState& state) = 0;
};

} // namespace FDM
} // namespace SF
