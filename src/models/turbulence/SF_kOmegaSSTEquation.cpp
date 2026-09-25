/// @file SF_kOmegaSSTEquation.cpp
/// @brief 可注册湍流输运方程与模型操作实现。

#include "SF_equationModels.h"

#include "SF_equationModelOps.h"
#include "SF_phaseProperties.h"
#include "methods/numerics/structured/SF_vectorCalculus.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace SF::Turbulence {
namespace {

struct SSTConstants {
    double sigmaK1=0.85;
    double sigmaK2=1.0;
    double sigmaW2=0.856;
    double beta2=0.0828;
    double gamma2=0.44;
};

} // namespace

void prepareKOmegaSST(
        const Physics::PhaseSystems::PhaseSystem& system,
        const FDM::TurbulenceConfig& config,
        PhaseEquationState& state) {
    const auto& c=config.coefficients;
    const SSTConstants fixed;
    if (!(std::isfinite(c.betaStar)&&c.betaStar>0.0
          &&std::isfinite(c.beta1)&&c.beta1>0.0
          &&std::isfinite(c.gamma1)&&c.gamma1>0.0
          &&std::isfinite(c.a1)&&c.a1>0.0
          &&std::isfinite(c.sigmaOmega1)&&c.sigmaOmega1>0.0)) {
        throw std::runtime_error(
            "kOmegaSST coefficients must be finite and positive.");
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
            ModelOps::clearEquationCell(state.specificDissipation,cell);
            continue;
        }
        const double alpha=phase.primitive.alpha.values()[(size_t)cell];
        const double rho=phase.primitive.density.values()[(size_t)cell];
        const double mass=phase.primary.phaseMass.values()[(size_t)cell];
        const double kValue=
            state.kineticEnergy.variable.values()[(size_t)cell];
        const double omega=
            state.specificDissipation.variable.values()[(size_t)cell];
        ModelOps::requirePositiveState(
            "kOmegaSST",phase.name,"omega",rho,mass,kValue,omega,
            c.kFloor,c.omegaFloor,i,j,k);
        const double distance=state.wallDistance.values()[(size_t)cell];
        const bool wallPoint=std::isfinite(distance)&&distance==0.0;
        const double temperature=
            phase.primitive.temperature.values()[(size_t)cell];
        const double mu=Physics::Multiphase::phaseViscosity(
            properties,temperature);
        if (wallPoint) {
            state.eddyViscosity.values()[(size_t)cell]=0.0;
            state.kineticEnergy.diffusivity.values()[(size_t)cell]=
                alpha*mu;
            state.specificDissipation.diffusivity.values()[(size_t)cell]=
                alpha*mu;
            state.kineticEnergy.explicitSource.values()[(size_t)cell]=0.0;
            state.kineticEnergy.implicitSink.values()[(size_t)cell]=0.0;
            state.specificDissipation.explicitSource
                .values()[(size_t)cell]=0.0;
            state.specificDissipation.implicitSink
                .values()[(size_t)cell]=0.0;
            continue;
        }
        if (!std::isfinite(distance)||distance<=0.0) {
            throw std::runtime_error(
                "kOmegaSST requires positive wall distance away from walls.");
        }
        const auto gradK=CENTRAL2::grad(
            geometry,state.kineticEnergy.variable,i,j,k);
        const auto gradW=CENTRAL2::grad(
            geometry,state.specificDissipation.variable,i,j,k);
        const double crossDot=SF::dot(gradK,gradW);
        // SST 的 CD_kω 下界属于模型定义，用于构造连续且有界的混合函数。
        const double crossCoefficient=std::max(
            2.0*rho*fixed.sigmaW2*crossDot/omega,1.0e-20);
        const double nu=mu/rho;
        const double rootK=std::sqrt(kValue);
        const double firstArgument=std::min(
            std::max(
                rootK/(c.betaStar*omega*distance),
                500.0*nu/(distance*distance*omega)),
            4.0*rho*fixed.sigmaW2*kValue
                /(crossCoefficient*distance*distance));
        const double f1=std::tanh(std::pow(firstArgument,4.0));
        const double secondArgument=std::max(
            2.0*rootK/(c.betaStar*omega*distance),
            500.0*nu/(distance*distance*omega));
        const double f2=std::tanh(
            std::pow(secondArgument,2.0));
        const double vorticity=ModelOps::vorticityMagnitude(
            geometry,phase,i,j,k);
        // F2 与涡量限制器共同限制近壁和强剪切区的湍流黏度。
        const double muT=rho*c.a1*kValue
            /std::max(c.a1*omega,vorticity*f2);
        if (!std::isfinite(muT)||muT<0.0) {
            throw std::runtime_error(
                "kOmegaSST produced invalid turbulent viscosity.");
        }
        const double sigmaK=
            f1*fixed.sigmaK1+(1.0-f1)*fixed.sigmaK2;
        const double sigmaW=
            f1*c.sigmaOmega1+(1.0-f1)*fixed.sigmaW2;
        const double beta=f1*c.beta1+(1.0-f1)*fixed.beta2;
        const double gamma=f1*c.gamma1+(1.0-f1)*fixed.gamma2;
        state.eddyViscosity.values()[(size_t)cell]=muT;
        state.kineticEnergy.diffusivity.values()[(size_t)cell]=
            alpha*(mu+sigmaK*muT);
        state.specificDissipation.diffusivity.values()[(size_t)cell]=
            alpha*(mu+sigmaW*muT);
        const double production=alpha*muT
            *ModelOps::strainSquared(geometry,phase,i,j,k);
        // 两个耗散项进入矩阵对角；生产与交叉扩散保留为显式源。
        const double limitedProduction=std::min(
            production,20.0*c.betaStar*mass*kValue*omega);
        state.kineticEnergy.explicitSource.values()[(size_t)cell]=
            limitedProduction;
        state.kineticEnergy.implicitSink.values()[(size_t)cell]=
            c.betaStar*mass*omega;
        const double crossSource=2.0*(1.0-f1)*mass
            *fixed.sigmaW2*crossDot/omega;
        state.specificDissipation.explicitSource
            .values()[(size_t)cell]=
                gamma*rho/muT*limitedProduction+crossSource;
        state.specificDissipation.implicitSink
            .values()[(size_t)cell]=beta*mass*omega;
    }
}

} // namespace SF::Turbulence
