/// @file SF_DNS.cpp
/// @brief DNS 无湍流闭式分支及能力声明。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.05.28-----------*/

#include "SF_DNS.h"

namespace SF {
namespace Turbulence {
namespace DNS {

void DirectNumericalSimulationModel::initialize(const Field& flow,
                                                ScalarFields& state,
                                                const FDM::TurbulenceConfig& config) {
    correct(flow, state, config, 0.0);
}

void DirectNumericalSimulationModel::applyBoundary(const Field&,
                                                   ScalarFields&,
                                                   const FDM::TurbulenceConfig&) {
}

void DirectNumericalSimulationModel::correct(const Field& flow,
                                             ScalarFields& state,
                                             const FDM::TurbulenceConfig&,
                                             double) {
    int ng = flow.NG();
    int nx = flow.NX();
    int ny = flow.NY();
    int nz = flow.NZ();

    for (int k = ng; k < nz + ng; ++k)
        for (int j = ng; j < ny + ng; ++j)
            for (int i = ng; i < nx + ng; ++i)
                state.EddyMu(i, j, k) = 0.0;
}

double DirectNumericalSimulationModel::eddyDynamicViscosity(const Field&,
                                                            const ScalarFields&,
                                                            const FDM::TurbulenceConfig&,
                                                            int, int, int) const {
    return 0.0;
}

} // namespace DNS
} // namespace Turbulence
} // namespace SF
