#pragma once
#include "solver/run/SF_planExecutor.h"
#include "solver/system/SF_stateRealizer.h"
#include "solver/system/SF_numericalSystem.h"
#include "core/state/SF_stateBundle.h"
namespace SF::Time {
/// @brief Stable numerical callbacks; the compiled plan owns phase ordering.
struct TemporalCallbacks {
    std::string identity;
    /// @brief 启动阶段：所有 provider 先绑定数组，再绑定跨 provider 输入；不执行数值。
    std::function<void()> bindReads;
    std::function<void()> snapshot, validatePublish, publish, publishCoefficients;
    std::function<double()> stepSize;
    std::function<void(const Run::ExecutionContext&)> prepare, coefficients, rhs, advance;
};
void bindStageGroup(Run::OpRegistry&,const System::CompiledSolvePlan&,
    const System::CompiledNumericalSystem&,State::StateBundle&,
    System::StateRealization&,const double& limit,std::vector<TemporalCallbacks>,
    std::vector<std::function<void()>> inputBindings={});
}
