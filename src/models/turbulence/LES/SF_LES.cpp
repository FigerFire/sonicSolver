/// @file SF_LES.cpp
/// @brief LES 湍流模型选择与推进实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.05.28-----------*/

#include "SF_LES.h"

#include <algorithm>

namespace SF {
namespace Turbulence {
namespace LES {

void SmagorinskyModel::initialize(const Field& flow,
                                  ScalarFields& state,
                                  const FDM::TurbulenceConfig& config) {
    correct(flow, state, config, 0.0);
}

void SmagorinskyModel::applyBoundary(const Field&,
                                     ScalarFields&,
                                     const FDM::TurbulenceConfig&) {
}

void SmagorinskyModel::correct(const Field& flow,
                               ScalarFields& state,
                               const FDM::TurbulenceConfig& config,
                               double) {
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

double SmagorinskyModel::eddyDynamicViscosity(const Field& flow,
                                              const ScalarFields&,
                                              const FDM::TurbulenceConfig& config,
                                              int i, int j, int k) const {
    double rho = std::max(flow(i, j, k, RHO), 1.0e-12);
    double delta = characteristicFilterWidth(flow, i, j, k, config.coefficients.filterScale);
    double sMag = strainRateMagnitude(flow, i, j, k);
    double cs = config.coefficients.cSmagorinsky;
    return rho * (cs * delta) * (cs * delta) * sMag;
}

} // namespace LES
} // namespace Turbulence
} // namespace SF
