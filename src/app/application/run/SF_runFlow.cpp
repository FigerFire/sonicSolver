/// @file SF_runFlow.cpp
/// @brief 统一的单场时间循环驱动器。

#include "SF_runFlow.h"

#include <stdexcept>
#include <utility>

namespace SF::Application::Runners {

int runFlow(
        FDM::INavierStokesStepper& stepper,
        FDM::SolverState& state,
        const Workflow::Plan& workflow,
        const Time::RunControl& control,
        const std::string& stepperName,
        std::function<void(int)> saveStep,
        std::function<void(double)> saveTime,
        std::function<void(int, const Time::AdvanceResult&)> report,
        std::function<bool(bool)> globallyFinished,
        std::function<std::string(const FDM::StepResult&)> formatDetail) {
    if (!formatDetail) {
        throw std::runtime_error("runFlow requires a detail formatter.");
    }
    if (workflow.stages.empty()) {
        throw std::runtime_error("runFlow requires a non-empty workflow plan.");
    }
    stepper.bindSolveStages(workflow.stages);
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

} // namespace SF::Application::Runners
