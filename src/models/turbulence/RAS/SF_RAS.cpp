/// @file SF_RAS.cpp
/// @brief RAS 湍流模型选择与闭式实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.05.28-----------*/

#include "SF_RAS.h"

namespace SF {
namespace Turbulence {
namespace RAS {

void RASModelBase::initializeScalar(const Field& flow,
                                    ScalarFields& state,
                                    const std::vector<BCSetting<double>>& settings,
                                    ScalarSlot slot) const {
    applyScalarInitialConditions(flow, settings, state, slot);
}

void RASModelBase::applyScalarBC(const Field& flow,
                                 ScalarFields& state,
                                 const std::vector<BCSetting<double>>& settings,
                                 ScalarSlot slot) const {
    applyScalarBoundaryConditions(flow, settings, state, slot);
}

void RASModelBase::refreshEddyMu(const Field& flow,
                                 ScalarFields& state,
                                 const FDM::TurbulenceConfig& config) const {
    int ng = flow.NG();
    int nx = flow.NX();
    int ny = flow.NY();
    int nz = flow.NZ();

    for (int k = ng; k < nz + ng; ++k) {
        for (int j = ng; j < ny + ng; ++j) {
            for (int i = ng; i < nx + ng; ++i) {
                if (flow.CellFlag(i, j, k) != FLUID_CELL) {
                    state.EddyMu(i, j, k) = 0.0;
                    continue;
                }
                state.EddyMu(i, j, k) = eddyDynamicViscosity(flow, state, config, i, j, k);
            }
        }
    }
}

} // namespace RAS
} // namespace Turbulence
} // namespace SF
