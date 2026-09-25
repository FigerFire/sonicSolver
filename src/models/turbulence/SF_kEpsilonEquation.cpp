/// @file SF_kEpsilonEquation.cpp
/// @brief 可注册湍流输运方程与模型操作实现。

#include "SF_equationModels.h"

#include "SF_equationModelOps.h"
#include "SF_phaseProperties.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace SF::Turbulence {

void prepareKEpsilon(
        const Physics::PhaseSystems::PhaseSystem& system,
        const FDM::TurbulenceConfig& config,
        PhaseEquationState& state) {
    const auto& c=config.coefficients;
    if (!(std::isfinite(c.cMu)&&c.cMu>0.0
          &&std::isfinite(c.c1)&&c.c1>0.0
          &&std::isfinite(c.c2)&&c.c2>0.0
          &&std::isfinite(c.sigmaK)&&c.sigmaK>0.0
          &&std::isfinite(c.sigmaEpsilon)&&c.sigmaEpsilon>0.0)) {
        throw std::runtime_error(
            "kEpsilon coefficients must be finite and positive.");
    }
    const Field& geometry=system.geometry();
    const auto& phase=system.phases()[state.phaseIndex];
    const auto& properties=system.phaseProperties(state.phaseIndex);
    for (int cell=0;cell<geometry.TotalSize();++cell) {
        int i=0,j=0,k=0;
        geometry.getIJK(cell,i,j,k);
        if (!ModelOps::isPhysical(geometry,i,j,k)) {
            state.eddyViscosity.values()[(size_t)cell]=0.0;
            ModelOps::clearEquationCell(state.kineticEnergy,cell);
            ModelOps::clearEquationCell(state.dissipation,cell);
            continue;
        }
        const double alpha=phase.primitive.alpha.values()[(size_t)cell];
        const double rho=phase.primitive.density.values()[(size_t)cell];
        const double mass=phase.primary.phaseMass.values()[(size_t)cell];
        const double kValue=
            state.kineticEnergy.variable.values()[(size_t)cell];
        const double epsilon=
            state.dissipation.variable.values()[(size_t)cell];
        ModelOps::requirePositiveState(
            "kEpsilon",phase.name,"epsilon",rho,mass,kValue,epsilon,
            c.kFloor,c.epsilonFloor,i,j,k);
        const double muT=c.cMu*rho*kValue*kValue/epsilon;
        if (!std::isfinite(muT)||muT<0.0) {
            throw std::runtime_error(
                "kEpsilon produced invalid turbulent viscosity.");
        }
        state.eddyViscosity.values()[(size_t)cell]=muT;
        const double temperature=
            phase.primitive.temperature.values()[(size_t)cell];
        const double mu=Physics::Multiphase::phaseViscosity(
            properties,temperature);
        state.kineticEnergy.diffusivity.values()[(size_t)cell]=
            alpha*(mu+muT/c.sigmaK);
        state.dissipation.diffusivity.values()[(size_t)cell]=
            alpha*(mu+muT/c.sigmaEpsilon);
        // 生产项显式装配，耗散项作为非负对角汇半隐式装配。
        const double production=alpha*muT
            *ModelOps::strainSquared(geometry,phase,i,j,k);
        const double limitedProduction=
            std::min(production,50.0*mass*epsilon);
        state.kineticEnergy.explicitSource.values()[(size_t)cell]=
            limitedProduction;
        state.kineticEnergy.implicitSink.values()[(size_t)cell]=
            mass*epsilon/kValue;
        state.dissipation.explicitSource.values()[(size_t)cell]=
            c.c1*epsilon/kValue*limitedProduction;
        state.dissipation.implicitSink.values()[(size_t)cell]=
            c.c2*mass*epsilon/kValue;
    }
}

} // namespace SF::Turbulence
