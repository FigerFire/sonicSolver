#pragma once
#include "solver/run/SF_planExecutor.h"
#include "solver/system/SF_stateRealizer.h"
#include "solver/system/SF_numericalSystem.h"
#include "core/state/SF_stateBundle.h"
namespace SF::Time {
/// @brief Stable numerical callbacks; the compiled plan owns phase ordering.
struct TemporalCallbacks {
    std::string identity;
    std::function<void()> snapshot, validatePublish, publish;
    std::function<double()> stepSize;
    std::function<void(const Run::ExecutionContext&)> prepare, rhs, advance;
};
void bindStageGroup(Run::OpRegistry&,const System::CompiledSolvePlan&,
    const System::CompiledNumericalSystem&,State::StateBundle&,
    System::StateRealization&,const double& limit,std::vector<TemporalCallbacks>);
}
