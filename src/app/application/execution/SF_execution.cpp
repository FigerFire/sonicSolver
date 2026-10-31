/// @file SF_execution.cpp
/// @brief 按已请求 STATE backing 和操作资源装配；数值顺序仅来自冻结计划。
/// 相与单流体的真实存储/发布需求分别装配，不作为缺失资源的 fallback。

#include "app/application/execution/SF_execution.h"
#include "app/application/SF_application.h"
#include "SF_environment.h"
#include "app/application/adapters/SF_scalarResources.h"
#include "app/application/output/SF_fields.h"
#include "core/state/SF_stateBundle.h"
#include "infrastructure/execution/SF_localExecutionRuntime.h"
#include "core/interfaces/SF_log.h"
#include "solver/system/SF_stateRealizer.h"
#include "solver/system/SF_systemValidator.h"

#include <stdexcept>

namespace SF::Application::Execution {

int executeSingle(
        Field& field, ResultWriter& writer,
        Parallel::ParallelContext& parallel, IBM::IB& ibm,
        const FDM::SolverConfig& solverConfig,
        const System::ResolvedSimulationSystem& system,
        const System::CompiledSolvePlan& plan,
        const CaseConfig& caseConfig, bool ibmEnabled,
        bool initialOutputOnly, int localBlockId) {
    if (&plan!=&system.solvePlan)
        throw std::runtime_error("Resource assembly requires the authoritative compiled plan identity.");
    System::validate(system);
    if (System::requiredPhaseInputs(System::programOf(system))) {
        return Detail::withPhaseResources(field,writer,parallel,solverConfig,system,plan,
            caseConfig,initialOutputOnly,localBlockId);
    }
    const auto run=[&](State::StateBundle& state,FDM::SolverServices services,
                       const Time::DriverCallbacks& callbacks) {
        return executeEquations(solverConfig,system,state,std::move(services),caseConfig.time,callbacks);
    };
    const int packed=System::requiredPackedComponents(system.executableSystem.state,plan.compiledProgram,"conservative");
    const auto* pressure=System::requiredPressureInputs(System::programOf(system));
    if (packed || pressure) {
        if (packed && field.NVar()<packed)
            throw std::runtime_error("Insufficient packed STATE backing: required="+std::to_string(packed)
                +" available="+std::to_string(field.NVar()));
        return Detail::withSingleFieldResources(
            field,writer,parallel,ibm,solverConfig,system,plan,caseConfig,
            ibmEnabled,initialOutputOnly,localBlockId,run);
    }
    State::StateBundle state;state.patches={&field};state.time=caseConfig.time.startTime;
    auto resources=Detail::bindScalarResources(field,state,system,caseConfig);
    if (resources.instances.empty() || field.NVar()!=0 || field.stateModel())
        throw std::runtime_error("Compiled system has no implemented physical-resource assembly for its frozen requirements.");
    for (const auto& symbol:system.executableSystem.state.symbols())
        if (!symbol.constantValue && !symbol.producer
            && state.distributed.select(symbol.storageKey,State::HaloSyncStage::None).empty())
            throw std::runtime_error("Unbound physical STATE: "+symbol.id);
    FDM::SolverServices services;::SF::Execution::LocalRuntime runtime;
    services.executionRuntime=&runtime;services.scalarInstances=resources.instances;
    const auto output=Output::registeredStateFields(system.executableSystem.state,state,field);
    auto saveStep=[&](int step){writer.saveStep(field,step,state.time,output);};
    auto saveTime=[&](double time){writer.save(field,time,output);};
    SF::broadcast("Native resources: ","geometry-only, Q=0, EOS=none, scalar arrays="+std::to_string(resources.arrays.size()));
    if(caseConfig.writeInitial || initialOutputOnly) {if(caseConfig.time.writeByStep)saveStep(0);else saveTime(state.time);}
    if(initialOutputOnly)return 0;
    Time::DriverCallbacks callbacks;callbacks.saveStep=saveStep;callbacks.saveTime=saveTime;
    callbacks.report=[&](int step,const Time::AdvanceResult& result){SF::broadcast("Scalar step: ","step="+std::to_string(step)+", time="+std::to_string(result.time)+", dt="+std::to_string(result.dt));};
    return run(state,std::move(services),callbacks);
}

int execute(const System::ResolvedSimulationSystem& system,
            const System::CompiledSolvePlan& plan,
            const CaseConfig& caseConfig,
            ExecutionEnvironment& env,
            const RunRequest& request) {
    if (env.domain == DomainKind::DistributedMultiPatch) {
        if (system.executableSystem.immersed
            && system.executableSystem.immersed->wallClosure==FDM::ImmersedWallClosure::StationaryNoSlipAdiabatic)
            throw std::runtime_error("stationaryNoSlipAdiabatic requires a single physical patch; composite/distributed wall ports are Unsupported.");
        return executeMulti(
            env.parallelMesh, *env.writer, *env.parallel, env.compositeIBM,
            caseConfig.solver, system, plan, caseConfig, env.ibmEnabled,
            request.initialOutputOnly);
    }
    return executeSingle(
        *env.activeField, *env.writer, *env.parallel, env.ibm,
        caseConfig.solver, system, plan, caseConfig, env.ibmEnabled,
        request.initialOutputOnly, env.localBlockId);
}

} // namespace SF::Application::Execution
