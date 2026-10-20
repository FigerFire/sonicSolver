/// @file SF_smagorinskyEquation.cpp
/// @brief 可注册湍流输运方程与模型操作实现。

#include "SF_equationModels.h"

#include "SF_equationModelOps.h"

#include <cmath>
#include <stdexcept>

namespace SF::Turbulence {

void prepareSmagorinsky(
        const Physics::PhaseSystems::PhaseSystem& system,
        const FDM::TurbulenceConfig& config,
        PhaseEquationState& state) {
    const double coefficient=config.coefficients.cSmagorinsky;
    if (!std::isfinite(coefficient)||coefficient<=0.0) {
        throw std::runtime_error(
            "Smagorinsky coefficient must be finite and positive.");
    }
    const Field& geometry=system.geometry();
    const auto& phase=system.phases()[state.phaseIndex];
    for (int cell=0;cell<geometry.TotalSize();++cell) {
        int i=0,j=0,k=0;
        geometry.getIJK(cell,i,j,k);
        if (!ModelOps::isPhysical(geometry,i,j,k)) {
            state.eddyViscosity.values()[(size_t)cell]=0.0;
            continue;
        }
        const double rho=
            phase.primitive.density.values()[(size_t)cell];
        if (!std::isfinite(rho)||rho<=0.0) {
            throw std::runtime_error(
                "Smagorinsky received invalid phase density.");
        }
        const double delta=ModelOps::filterWidth(
            geometry,i,j,k,config.coefficients.filterScale);
        const double strain=std::sqrt(
            ModelOps::strainSquared(geometry,phase,i,j,k));
        // LES 只更新代数 mu_t，不创建伪造的 k/epsilon 输运状态。
        const double muT=rho
            *(coefficient*delta)*(coefficient*delta)*strain;
        if (!std::isfinite(muT)||muT<0.0) {
            throw std::runtime_error(
                "Smagorinsky produced invalid turbulent viscosity.");
        }
        state.eddyViscosity.values()[(size_t)cell]=muT;
    }
}

} // namespace SF::Turbulence
