/// @file test_explicitStages.cpp
/// @brief Freezes SSP-RK3 stage times and state combinations after Plan lowering.

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

int main() {
    SF::Field field;
    field.setup(1,1,1,0,5);
    for (int variable = 0; variable < field.NVar(); ++variable) {
        field(0,0,0,variable) = 10.0;
    }

    std::vector<SF::Field*> fields{&field};
    std::vector<SF::SolverAlgorithm::PatchWorkspace> workspaces(1);
    workspaces.front().ensureFor(field);
    SF::State::StateBundle state;
    state.patches = fields;
    state.time = 2.0;
    state.dt = 0.1;

    std::vector<double> stageTimes;
    int publishCount = 0;
    int validateCount = 0;
    SF::Time::Explicit::Workspace workspace;
    SF::Time::Explicit::begin(
        workspace,fields,state,SF::FDM::TimeScheme::SSPRK3,nullptr);
    for (int stage = 0; stage < 3; ++stage) {
        SF::Time::Explicit::executeStage(
            workspace,stage,fields,workspaces,state,nullptr,
            [&](const std::vector<SF::Field*>&,
                std::vector<SF::SolverAlgorithm::PatchWorkspace>& rhs,
                double stageTime) {
                stageTimes.push_back(stageTime);
                rhs.front().residual.clear();
                for (int variable = 0; variable < field.NVar(); ++variable) {
                    rhs.front().residual.source(0,0,0,variable) = -2.0;
                }
            },
            [&](const std::vector<SF::Field*>&) { ++publishCount; },
            [&](const std::vector<SF::Field*>&, const char*) {
                ++validateCount;
            });
    }

    require(stageTimes.size() == 3,"wrong SSP-RK3 RHS count");
    require(std::abs(stageTimes[0]-2.0) < 1e-14
            && std::abs(stageTimes[1]-2.1) < 1e-14
            && std::abs(stageTimes[2]-2.05) < 1e-14,
            "SSP-RK3 stage times changed");
    require(std::abs(field(0,0,0,0)-9.8) < 1e-14,
            "SSP-RK3 Shu-Osher combination changed");
    require(publishCount == 3 && validateCount == 3,
            "SSP-RK3 publish/validate positions changed");
    require(!workspace.active && workspace.nextStage == 3,
            "explicit workspace did not finish after the final stage");
    return 0;
}
