#pragma once

/// @file SF_flowLoop.h
/// @brief 统一的单场时间循环驱动器（Generic FlowRunner 的串行/单场入口）。

#include "core/interfaces/SF_solverStepper.h"
#include "solver/system/SF_resolvedSimulationSystem.h"
#include "solver/algorithm/time/SF_time.h"

#include <functional>
#include <string>

namespace SF::Application::Execution {

/// @brief 运行一次完整的单场时间推进。
///
/// 只拥有“时间循环”这一件事：把 RunControl 约束下的 step 上限、输出对齐、
/// 跨 rank 终止归约交给 Time::Driver，并把每个 stepper 返回的 StepResult
/// 翻译成 Driver 需要的 AdvanceResult。
///
/// 方程装配、边界、状态与 MPI 语义全部由调用方已经完成的 setup 与传入的
/// 回调负责；本函数不 bindServices、不做初始输出，也不触碰 physical state。
///
/// @param stepper        已 bindServices 的 NS 推进器。
/// @param state          求解状态（其 bundle 已挂好）。
/// @param plan           已解析的唯一结构化 solve plan；在时间循环前稳定绑定。
/// @param control        时间控制（startTime/endTime/endStep/write*）。
/// @param stepperName    错误信息里使用的名称。
/// @param saveStep       按 step 输出回调（可为空）。
/// @param saveTime       按 time 输出回调（可为空）。
/// @param report         每步报告回调（可为空）。
/// @param globallyFinished 跨 rank 终止归约（可为空，空则本地判定）。
/// @param formatDetail   把 StepResult 翻译成 AdvanceResult.detail。
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
    std::function<std::string(const FDM::StepResult&)> formatDetail);

namespace Detail {

/// Compose the shared conservative solver/state/service tail used by both
/// single- and multi-patch application setup. Topology-specific setup remains
/// with the caller; numerical lifecycle authority remains in the compiled plan.
int runConservative(
    const FDM::SolverConfig& config,
    const System::ResolvedSimulationSystem& system,
    const System::CompiledSolvePlan& plan,
    State::StateBundle& state,
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
    std::function<bool(bool)> globallyFinished);

} // namespace Detail

} // namespace SF::Application::Execution
