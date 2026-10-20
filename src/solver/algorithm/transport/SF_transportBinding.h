#pragma once
#include "core/interfaces/SF_solverStepper.h"
#include "solver/run/SF_planExecutor.h"
#include "solver/system/SF_runtimeRequirements.h"
namespace SF::SolverAlgorithm {
/// @brief 初始化既有输运数组的通用边界及闭合；不分配physical state或推进时间。
void initializeTransport(const Field&,FDM::ITransportModel&);
/// @brief 通用标量BC后执行模型专属floor/wall关系。
void applyTransportBoundary(const Field&,FDM::ITransportModel&);
/// @brief 原位once-step数值sweep与闭合发布；只能由已编译操作调度。
void updateOnceStepTransport(const Field&,FDM::ITransportModel&,double dt);
void bindTransportOperations(Run::OpRegistry&, const System::RuntimeRequirements&,
                             State::StateBundle&, const FDM::SolverServices&);
}
