#pragma once
#include "SF_scalarBinding.h"
namespace SF::SolverAlgorithm {
/// @brief 绑定冻结温度方程的数值操作；不拥有物理时间或第二份物理数组。
std::vector<Time::TemporalCallbacks> thermalParticipants(const System::CompiledNumericalSystem&,
    const System::CompiledSolvePlan&,State::StateBundle&,const FDM::SolverServices&,
    System::StateRealization&);
}
