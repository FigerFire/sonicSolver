/// @file SF_flowLoop.cpp
/// @brief 统一的单场时间循环驱动器。

#include "app/application/execution/SF_flowLoop.h"

#include "solver/algorithm/SF_singleFluidStepper.h"

#include <stdexcept>
#include <utility>

// CompiledSolvePlan
//       │
//       ▼
// stepper.bindSolvePlan(plan)
//       │
//       ▼
// stepper.prepare(state)
//       │
//       ▼
// Time::Driver::run(...)
//       │
//       ├── 给 state 设置本步最大 dt
//       │
//       ├── stepper.advance(state)
//       │
//       ├── 检查 step 是否接受
//       │
//       └── 返回 dt / time / detail
//       │
//       ▼
// save / report / finished

namespace SF::Application::Execution {

int flowLoop(
        FDM::INavierStokesStepper& stepper,
        FDM::SolverState& state,
        const System::CompiledSolvePlan& plan,
        const Time::RunControl& control,
        const std::string& stepperName,
        std::function<void(int)> saveStep,
        std::function<void(double)> saveTime,
        std::function<void(int, const Time::AdvanceResult&)> report,
        std::function<bool(bool)> globallyFinished,
        std::function<std::string(const FDM::StepResult&)> formatDetail) {
    if (!formatDetail) {
        throw std::runtime_error("flowLoop requires a detail formatter.");
    }
    if (plan.root.children.empty()) {
        throw std::runtime_error("flowLoop requires a non-empty compiled solve plan.");
    }

    // 把已经编译好的 mathematical solve plan 给 stepper。
    stepper.bindSolvePlan(plan);
    // Resolve all mathematical unknowns against stable storage before the
    // first timestep so missing state/workspace is an initialization error.
    stepper.prepare(state);

    Time::Driver(control).run({
        [&](double limit, int currentStep) {
            (void)currentStep;
            state.maximumTimeStep = limit;
            const FDM::StepResult result = stepper.advance(state);
            if (!result.accepted) {
                throw std::runtime_error(
                    stepperName + " stepper rejected a step: "
                    + result.message);
            }
            return Time::AdvanceResult{
                result.dt, result.time, formatDetail(result)};
        },
        std::move(saveStep),
        std::move(saveTime),
        std::move(report),
        std::move(globallyFinished)});

    return 0;
}

} // namespace SF::Application::Execution

namespace SF::Application::Execution::Detail {

int runConservative(
        const FDM::SolverConfig& config,
        const System::ResolvedSimulationSystem& system,
        const System::CompiledSolvePlan& plan,
        State::StateBundle& bundle,
        FDM::IBoundaryPipeline& boundary,
        FDM::IExecutionRuntime& runtime,
        FDM::ISolverObserver& observer,
        FDM::ImmersedCouplingPorts immersed,
        FDM::ITransportModel* transport,
        FDM::IEquationSystemCoupling* equations,
        const Time::RunControl& control,
        const std::string& name,
        std::function<void(int)> saveStep,
        std::function<void(double)> saveTime,
        std::function<bool(bool)> globallyFinished) {
    SolverAlgorithm::SingleFluidStepper solver(
        config,system.executableSystem,system.numericalSystem,
        system.solvePlan,system.runtime);
    FDM::SolverState state;
    state.bundle = &bundle;
    FDM::SolverServices services;
    services.boundaryPipeline = &boundary;
    services.executionRuntime = &runtime;
    services.observer = &observer;
    services.immersed = immersed;
    services.transportModel = transport;
    services.equationSystem = equations;
    solver.bindServices(services);
    return flowLoop(
        solver,state,plan,control,name,
        std::move(saveStep),std::move(saveTime),{},
        std::move(globallyFinished),
        [](const FDM::StepResult& result) { return result.message; });
}

} // namespace SF::Application::Execution::Detail

