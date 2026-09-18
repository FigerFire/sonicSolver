/// @file test_explicitStages.cpp
/// @brief Freezes Euler/SSP-RK3/RK4 stage mathematics after Plan lowering.

#include "solver/algorithm/time/SF_explicit.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace {

void require(bool condition, const char* message) {
    if (condition) return;
    std::cerr << "Explicit stage test failed: " << message << '\n';
    std::exit(1);
}

} // namespace

struct Outcome {
    std::vector<double> stageTimes;
    double finalValue = 0.0;
    int publishCount = 0;
    int validateCount = 0;
    int nextStage = 0;
    bool active = false;
};

Outcome runConstantRhs(SF::FDM::TimeScheme scheme) {
    SF::Field field;
    field.setup(1,1,1,0,5);
    field(0,0,0,0) = 1.0;
    field(0,0,0,1) = 0.0;
    field(0,0,0,2) = 0.0;
    field(0,0,0,3) = 0.0;
    field(0,0,0,4) = 10.0;

    std::vector<SF::Field*> fields{&field};
    std::vector<SF::SolverAlgorithm::PatchWorkspace> workspaces(1);
    workspaces.front().ensureFor(field);
    SF::State::StateBundle state;
    state.patches = fields;
    state.time = 2.0;
    state.dt = 0.1;

    Outcome outcome;
    SF::Time::Explicit::Workspace workspace;
    SF::Time::Explicit::begin(
        workspace,fields,state,scheme,nullptr);
    for (int stage = 0; stage < SF::Time::Explicit::stageCount(scheme);
         ++stage) {
        SF::Time::Explicit::executeStage(
            workspace,stage,fields,workspaces,state,nullptr,
            [&](const std::vector<SF::Field*>&,
                std::vector<SF::SolverAlgorithm::PatchWorkspace>& rhs,
                double stageTime) {
                outcome.stageTimes.push_back(stageTime);
                rhs.front().residual.clear();
                rhs.front().residual.source(0,0,0,4) = -2.0;
            },
            [&](const std::vector<SF::Field*>&) { ++outcome.publishCount; },
            [&](const std::vector<SF::Field*>&, const char*) {
                ++outcome.validateCount;
            });
    }
    outcome.finalValue = field(0,0,0,4);
    outcome.nextStage = workspace.nextStage;
    outcome.active = workspace.active;
    return outcome;
}

int main() {
    const Outcome euler = runConstantRhs(SF::FDM::TimeScheme::Euler);
    require(euler.stageTimes.size() == 1
            && std::abs(euler.stageTimes[0]-2.0) < 1e-14,
            "Euler stage time changed");
    require(std::abs(euler.finalValue-9.8) < 1e-14,
            "Euler update changed");
    require(euler.publishCount == 1 && euler.validateCount == 1,
            "Euler publish/validate positions changed");

    const Outcome ssp = runConstantRhs(SF::FDM::TimeScheme::SSPRK3);

    require(ssp.stageTimes.size() == 3,"wrong SSP-RK3 RHS count");
    require(std::abs(ssp.stageTimes[0]-2.0) < 1e-14
            && std::abs(ssp.stageTimes[1]-2.1) < 1e-14
            && std::abs(ssp.stageTimes[2]-2.05) < 1e-14,
            "SSP-RK3 stage times changed");
    require(std::abs(ssp.finalValue-9.8) < 1e-14,
            "SSP-RK3 Shu-Osher combination changed");
    require(ssp.publishCount == 3 && ssp.validateCount == 3,
            "SSP-RK3 publish/validate positions changed");
    require(!ssp.active && ssp.nextStage == 3,
            "explicit workspace did not finish after the final stage");

    const Outcome rk4 = runConstantRhs(SF::FDM::TimeScheme::RK4);
    require(rk4.stageTimes.size() == 4,"wrong RK4 RHS count");
    require(std::abs(rk4.stageTimes[0]-2.0) < 1e-14
            && std::abs(rk4.stageTimes[1]-2.05) < 1e-14
            && std::abs(rk4.stageTimes[2]-2.05) < 1e-14
            && std::abs(rk4.stageTimes[3]-2.1) < 1e-14,
            "RK4 stage times changed");
    require(std::abs(rk4.finalValue-9.8) < 1e-14,
            "RK4 combination changed");
    require(rk4.publishCount == 4 && rk4.validateCount == 5,
            "RK4 publish/validate positions changed");
    return 0;
}
