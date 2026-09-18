/// @file SF_eulerian.cpp
/// @brief 组装 Eulerian phase state、pressure providers 与 compiled-plan execution。
///
/// Data flow:
///   Field + executable equations + phase configuration
///       -> phase state/workspace/providers
///       -> PressureStepper callbacks executed by runFlow
///
/// 本文件不拥有 PIMPLE loop、不修改 pressure 数学，也不选择 MPI semantics。

#include "app/application/execution/SF_runners.h"
#include "app/application/execution/SF_executionAssemblers.h"
#include "app/application/run/SF_runFlow.h"

#include "SF_resultWriter.h"
#include "app/application/output/SF_report.h"
#include "core/interfaces/SF_log.h"
#include "app/application/output/SF_fields.h"
#include "SF_phaseSystem.h"
#include "SF_pressureStepper.h"
#include "app/application/adapters/SF_ibmAdapters.h"
#include "solver/algorithm/time/SF_time.h"
#include "SF_parallelContext.h"
#include "core/state/SF_state.h"
#include "infrastructure/execution/SF_executionRuntime.h"

#include <cmath>
#include <stdexcept>
#include <vector>

namespace SF::Application::Runners::Detail {

int executeEulerianEquations(
        Field& field,
        ResultWriter& writer,
        Parallel::ParallelContext& parallel,
        const FDM::SolverConfig& solverConfig,
        const System::ResolvedSimulationSystem& resolvedSystem,
        const System::CompiledSolvePlan& plan,
        const CaseConfig& caseConfig,
        bool initialOutputOnly,
        int localBlockId) {
    using Report::broadcastSolverConfig;
    using Report::formatRunControl;
    using Report::formatTimeStepStatus;
    using Report::formatTimeValue;
    using Output::eulerianTurbulenceVTKScalars;
    using Output::eulerianVTKScalars;
    auto& coordinator = parallel.coordinator();
    Execution::Runtime executionRuntime(
        parallel.active() && parallel.size() > 1 ? &coordinator : nullptr);

    if (std::abs(caseConfig.multiPhase.initialPressure
                 - solverConfig.pressure.workflow.referencePressure)
        > 1.0e-10 * solverConfig.pressure.workflow.referencePressure) {
        broadcast("Fatal: ",
            "phaseProperties initialPressure must equal solverProperties "
            "referencePressure for the shared-pressure reference cell.");
        return -1;
    }

    auto phaseSystem = Physics::PhaseSystems::makePhaseSystem(
        caseConfig.multiPhase);
    phaseSystem->initialize(field);
    EulerianEulerian::PressureStepper pressureStepper(
        *phaseSystem, solverConfig, resolvedSystem);

    broadcast("PhaseSystem       : ",
              "EulerianEulerian (independent U/rho/T, shared p, N-1 alpha)");
    broadcastSolverConfig(solverConfig, true);
    broadcast("Linear algebra    : ",
              "HYPRE IJ/ParCSR + FlexGMRES/PCG + BoomerAMG");
    if (pressureStepper.turbulence().active()) {
        broadcast("Turbulence equations: ",
                  pressureStepper.turbulence().description());
    }
    broadcast("Run control: ", formatRunControl(caseConfig));

    auto scalars = [&]() {
        auto values = eulerianVTKScalars(*phaseSystem);
        auto turbulence = eulerianTurbulenceVTKScalars(
            pressureStepper.turbulence());
        values.insert(values.end(), turbulence.begin(), turbulence.end());
        return values;
    };
    auto saveStep = [&](int outputStep) {
        auto values = scalars();
        if (parallel.active()) writer.save(field, outputStep, localBlockId, values);
        else writer.save(field, outputStep, values);
    };
    auto saveTime = [&](double outputTime) {
        auto values = scalars();
        if (parallel.active()) writer.save(field, outputTime, localBlockId, values);
        else writer.save(field, outputTime, values);
    };
    if (initialOutputOnly) {
        if (caseConfig.time.writeByStep) saveStep(0);
        else saveTime(caseConfig.time.startTime);
        return 0;
    }
    if (caseConfig.writeInitial) {
        if (caseConfig.time.writeByStep) saveStep(0);
        else saveTime(caseConfig.time.startTime);
    }

    State::StateBundle stateBundle;
    stateBundle.patches = {&field};
    stateBundle.time = caseConfig.time.startTime;
    stateBundle.registerConservativeState();
    pressureStepper.registerState(stateBundle);
    FDM::SolverState solverState;
    solverState.bundle = &stateBundle;
    FDM::SolverServices solverServices;
    solverServices.executionRuntime = &executionRuntime;
    pressureStepper.bindServices(solverServices);

    // EulerianEulerian 的 detail 需要跨 rank 归约 summary；归约顺序与原来
    // 完全一致（advance -> accepted 检查 -> 归约 -> 格式化）。
    auto formatDetail = [&](const FDM::StepResult&) {
        auto summary = pressureStepper.lastSummary();
        if (parallel.active()) {
            summary.massImbalance = executionRuntime.globalSum(
                summary.massImbalance);
            summary.maxAlphaSumError = executionRuntime.globalMaximum(
                summary.maxAlphaSumError);
            summary.totalPhaseMass = executionRuntime.globalSum(
                summary.totalPhaseMass);
            summary.totalPhaseEnthalpy = executionRuntime.globalSum(
                summary.totalPhaseEnthalpy);
            summary.totalWallHeat = executionRuntime.globalSum(
                summary.totalWallHeat);
            summary.pressureResidual = executionRuntime.globalMaximum(
                summary.pressureResidual);
            summary.turbulenceResidual = executionRuntime.globalMaximum(
                summary.turbulenceResidual);
        }
        return "outer=" + std::to_string(summary.outerCorrectors)
            + ", pCorr=" + std::to_string(summary.pressureCorrections)
            + ", pIter=" + std::to_string(summary.pressureIterations)
            + ", pResidual=" + formatTimeValue(summary.pressureResidual)
            + ", maxAlphaSumError=" + formatTimeValue(summary.maxAlphaSumError)
            + ", phaseMass=" + formatTimeValue(summary.totalPhaseMass)
            + ", phaseEnthalpy=" + formatTimeValue(summary.totalPhaseEnthalpy)
            + ", wallHeat=" + formatTimeValue(summary.totalWallHeat)
            + ", HYPRE(rebuilds/solves)="
            + std::to_string(summary.pressureStructureRebuilds) + "/"
            + std::to_string(summary.pressureLinearSolves)
            + ", turbulenceIter="
            + std::to_string(summary.turbulenceIterations)
            + ", turbulenceResidual="
            + formatTimeValue(summary.turbulenceResidual)
            + ", turbulenceHYPRE(rebuilds/solves)="
            + std::to_string(summary.turbulenceStructureRebuilds) + "/"
            + std::to_string(summary.turbulenceLinearSolves);
    };
    auto report = [&](int, const Time::AdvanceResult& result) {
        broadcast("Eulerian step: ",
            formatTimeStepStatus(result.time, result.dt)
            + ", " + result.detail);
    };
    return runFlow(
        pressureStepper, solverState, plan, caseConfig.time,
        "Eulerian unified",
        saveStep,
        saveTime,
        report,
        [&](bool finished) { return coordinator.allRanksAgree(finished); },
        formatDetail);
}

} // namespace SF::Application::Runners::Detail
