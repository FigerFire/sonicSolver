/// @file SF_eulerian.cpp
/// @brief 应用持有相/模型/边界/输出；独立绑定函数借用资源执行原冻结计划。
///
/// Data flow:
///   Field + executable equations + phase configuration
///       -> phase state/workspace/providers
///       -> EulerianStepper callbacks executed by flowLoop
///
/// 本文件不拥有 PIMPLE loop、不修改 pressure 数学，也不选择 MPI semantics。

#include "app/application/execution/SF_execution.h"
#include "app/application/execution/SF_flowLoop.h"

#include "SF_resultWriter.h"
#include "app/application/output/SF_report.h"
#include "core/interfaces/SF_log.h"
#include "app/application/output/SF_fields.h"
#include "SF_phaseSystem.h"
#include "SF_eulerianStepper.h"
#include "app/application/adapters/SF_ibmAdapters.h"
#include "solver/algorithm/time/SF_time.h"
#include "SF_parallelContext.h"
#include "core/state/SF_state.h"
#include "infrastructure/execution/SF_executionRuntime.h"

#include "SF_phaseGravity.h"
#include "SF_phaseMRF.h"
#include "SF_phaseWallHeat.h"
#include "SF_phaseChange.h"
#include "solver/system/SF_eulerianAssembly.h"
#include "solver/system/SF_eulerianTurbulence.h"
#include <cmath>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <array>

namespace SF::Application::Execution::Detail {

namespace {
int runPhasePlan(Physics::PhaseSystems::PhaseSystem& phaseSystem,Turbulence::EquationSystem& turbulence,
        Physics::PhaseSystems::PhaseSourceRegistry& sources,Physics::PhaseSystems::PhaseBoundaryApplicator& boundary,
        State::StateBundle& stateBundle,FDM::SolverServices services,const System::ResolvedSimulationSystem& resolvedSystem,
        Parallel::ParallelContext& parallel,const Time::RunControl& control,const Time::DriverCallbacks& callbacks) {
    using Report::formatTimeValue;using Report::formatTimeStepStatus;
    auto& executionRuntime=*services.executionRuntime;
    EulerianEulerian::EulerianStepper pressureStepper(phaseSystem,turbulence,sources,boundary,
        resolvedSystem.executableSystem,resolvedSystem.numericalSystem,resolvedSystem.solvePlan,resolvedSystem.runtime);
    pressureStepper.registerState(stateBundle);
    pressureStepper.bindServices(services);
    FDM::SolverState solverState;solverState.bundle=&stateBundle;
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
    return flowLoop(
        pressureStepper, solverState, resolvedSystem.solvePlan, control,
        "Eulerian unified",
        callbacks.saveStep,
        callbacks.saveTime,
        report,
        [&](bool finished) { return parallel.coordinator().allRanksAgree(finished); },
        formatDetail);
}
}

int withPhaseResources(
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
    ::SF::Execution::Runtime executionRuntime(
        parallel.active() && parallel.size() > 1 ? &coordinator : nullptr);

    if(&plan!=&resolvedSystem.solvePlan)throw std::runtime_error("Phase resources require the authoritative plan identity.");
    const auto* inputs=System::requiredPhaseInputs(System::programOf(resolvedSystem));
    if(!inputs || resolvedSystem.runtime.report.status!=System::RuntimeStatus::Runnable)
        throw std::runtime_error("Phase resource assembly unsupported: "+resolvedSystem.runtime.report.reason);
    if (std::abs(caseConfig.multiPhase.initialPressure
                 - inputs->reference.referencePressure)
        > 1.0e-10 * inputs->reference.referencePressure) {
        broadcast("Fatal: ",
            "phaseProperties initialPressure must equal registered pressure "
            "referencePressure for the shared-pressure reference cell.");
        return -1;
    }

    auto phaseSystem = Physics::PhaseSystems::makePhaseSystem(
        caseConfig.multiPhase);
    phaseSystem->initialize(field);
    std::vector<std::string> phaseNames;
    for(const auto& phase:phaseSystem->phases())phaseNames.push_back(phase.name);
    System::validateEulerianAssemblyBindings(plan.compiledProgram,phaseNames,phaseSystem->referencePhaseIndex());
    Physics::PhaseSystems::PhaseBoundaryApplicator boundary;
    Physics::PhaseSystems::PhaseSourceRegistry sources;
    struct Source {
        FDM::SourceKind kind;System::EulerianSourceExtension extension;
        std::unique_ptr<Physics::PhaseSystems::PhaseEquationSource> (*make)(const FDM::SourceConfig&);
    };
    const std::array<Source,3> implementations{{
        {FDM::SourceKind::Gravity,System::EulerianSourceExtension::Gravity,[](const auto& c){return Physics::PhaseSystems::makePhaseGravitySource(c.gravity);}},
        {FDM::SourceKind::MRF,System::EulerianSourceExtension::MRF,[](const auto& c){return Physics::PhaseSystems::makePhaseMRFSource(c.rotating);}},
        {FDM::SourceKind::WallHeat,System::EulerianSourceExtension::WallHeat,[](const auto& c){return Physics::PhaseSystems::makePhaseWallHeatSource(c.wallHeat);}}
    }};
    for(const auto& call:plan.compiledProgram.steps)if(const auto* contract=System::eulerianAssemblyContract(call)) {
        for(const auto& source:implementations) {
            const bool consumer=source.kind==FDM::SourceKind::WallHeat
                ? contract->relation==System::EulerianAssemblyRelation::PhaseEnthalpyBalance
                : contract->relation==System::EulerianAssemblyRelation::PhaseMomentumBalance;
            if(!consumer)continue;
            const bool enabled=std::find(solverConfig.sources.enabled.begin(),solverConfig.sources.enabled.end(),source.kind)!=solverConfig.sources.enabled.end();
            const bool declared=std::find(contract->extensions.begin(),contract->extensions.end(),source.extension)!=contract->extensions.end();
            if(enabled!=declared)throw std::runtime_error("Phase source resource disagrees with mathematics: "+call.source.equation+" / "+System::toString(source.extension));
        }
    }
    bool wallHeatEnabled=false;
    for(const auto kind:solverConfig.sources.enabled) {
        const auto found=std::find_if(implementations.begin(),implementations.end(),[&](const auto& source){return source.kind==kind;});
        if(found==implementations.end())throw std::runtime_error("No phase source resource for selected SourceKind.");
        sources.add(found->make(solverConfig.sources));
        wallHeatEnabled=wallHeatEnabled || kind==FDM::SourceKind::WallHeat;
    }
    if(phaseSystem->config().phaseChange.enabled
        && Physics::PhaseChange::normalizeModel(phaseSystem->config().phaseChange.model)=="rpi" && !wallHeatEnabled)
        throw std::runtime_error("Eulerian RPI requires a wallHeatFlux boundary so the unified wall-boiling ledger can assemble mass, momentum and enthalpy.");
    Turbulence::EquationSystem turbulence(solverConfig.turbulence);
    turbulence.initialize(*phaseSystem);
    System::validateEulerianTurbulenceBindings(plan.compiledProgram,phaseNames,turbulence.config());

    broadcast("PhaseSystem       : ",
              "EulerianEulerian (independent U/rho/T, shared p, N-1 alpha)");
    broadcastSolverConfig(solverConfig, true);
    broadcast("Linear algebra    : ",
              "HYPRE IJ/ParCSR + FlexGMRES/PCG + BoomerAMG");
    if (turbulence.active()) {
        broadcast("Turbulence equations: ",
                  turbulence.description());
    }
    broadcast("Run control: ", formatRunControl(caseConfig));

    auto scalars = [&]() {
        auto values = eulerianVTKScalars(*phaseSystem);
        auto modelFields=eulerianTurbulenceVTKScalars(turbulence);
        values.insert(values.end(),modelFields.begin(),modelFields.end());
        return values;
    };
    State::StateBundle stateBundle;
    auto saveStep = [&](int outputStep) {
        auto values = scalars();
        writer.saveStep(field,outputStep,stateBundle.time,values,parallel.active()?localBlockId:-1);
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

    stateBundle.patches = {&field};
    stateBundle.time = caseConfig.time.startTime;
    stateBundle.registerConservativeState();
    FDM::SolverServices solverServices;solverServices.executionRuntime=&executionRuntime;
    Time::DriverCallbacks callbacks;callbacks.saveStep=saveStep;callbacks.saveTime=saveTime;
    return runPhasePlan(*phaseSystem,turbulence,sources,boundary,stateBundle,solverServices,
        resolvedSystem,parallel,caseConfig.time,callbacks);
}

} // namespace SF::Application::Execution::Detail
