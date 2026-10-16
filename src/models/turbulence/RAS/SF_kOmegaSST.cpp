/// @file SF_kOmegaSST.cpp
/// @brief RAS 湍流模型选择与闭式实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.05.28-----------*/

#include "SF_kOmegaSST.h"

#include "methods/numerics/structured/SF_structured.h"
#include "SF_viscous.h"
#include "methods/numerics/structured/SF_vectorCalculus.h"

#include <algorithm>
#include <cmath>

namespace SF {
namespace Turbulence {
namespace RAS {

namespace {

struct SSTConstants {
    double betaStar = 0.09;
    double sigmaK1 = 0.85;
    double sigmaW1 = 0.5;
    double beta1 = 0.075;
    double gamma1 = 5.0 / 9.0;
    double sigmaK2 = 1.0;
    double sigmaW2 = 0.856;
    double beta2 = 0.0828;
    double gamma2 = 0.44;
    double a1 = 0.31;
};

Vector3 scalarPhysicalGradient(const Field& flow,
                               const ScalarFields& state,
                               ScalarSlot slot,
                               int i, int j, int k) {
    return CENTRAL2::gradSampled(
        flow,
        [&](int ii, int jj, int kk) {
            return state.clamped(slot, ii, jj, kk);
        },
        i, j, k);
}

double vorticityMagnitude(const Field& flow, int i, int j, int k) {
    return norm(Math::curl(Viscous::velocityGradientAt(
        flow, i, j, k, Viscous::CentralOrder::SECOND)));
}

double laminarKinematicViscosity(double rho, double laminarMu) {
    return std::max(laminarMu, 0.0) / std::max(rho, 1.0e-12);
}

double wallDistance(const Field& flow, int i, int j, int k,bool strict) {
    const double d = flow.wallDistance(i, j, k);
    if (std::isfinite(d) && d > 1.0e-14) return d;
    if (strict) throw std::runtime_error("SST stationary immersed wall requires positive physical wall distance at ("
        +std::to_string(i)+","+std::to_string(j)+","+std::to_string(k)+"), d="+std::to_string(d)
        +", flag="+std::to_string(flow.CellFlag(i,j,k))+"; no far-field fallback.");
    return 1.0e20;
}

double crossDiffusionCD(const Field& flow,
                        const ScalarFields& state,
                        int i, int j, int k,
                        double rho,
                        double omega,
                        const SSTConstants& c) {
    const Vector3 gradK = scalarPhysicalGradient(
        flow, state, ScalarSlot::K, i, j, k);
    const Vector3 gradW = scalarPhysicalGradient(
        flow, state, ScalarSlot::Omega, i, j, k);
    return std::max(2.0 * rho * c.sigmaW2 * SF::dot(gradK, gradW)
                    / std::max(omega, 1.0e-30),
                    1.0e-20);
}

double blendingF1(const Field& flow,
                  const ScalarFields& state,
                  int i, int j, int k,
                  double rho,
                  double kVal,
                  double omega,
                  double laminarMu,
                  const SSTConstants& c,bool strict) {
    const double d = wallDistance(flow, i, j, k,strict);
    const double nu = laminarKinematicViscosity(rho, laminarMu);
    const double sqrtK = std::sqrt(std::max(kVal, 0.0));
    const double cd = crossDiffusionCD(flow, state, i, j, k, rho, omega, c);

    const double term1 = sqrtK / std::max(c.betaStar * omega * d, 1.0e-30);
    const double term2 = 500.0 * nu / std::max(d * d * omega, 1.0e-30);
    const double term3 = 4.0 * rho * c.sigmaW2 * kVal / std::max(cd * d * d, 1.0e-30);
    const double arg1 = std::min(std::max(term1, term2), term3);
    return std::tanh(std::pow(std::max(arg1, 0.0), 4.0));
}

double blendingF2(const Field& flow,
                  int i, int j, int k,
                  double rho,
                  double kVal,
                  double omega,
                  double laminarMu,
                  const SSTConstants& c,bool strict) {
    const double d = wallDistance(flow, i, j, k,strict);
    const double nu = laminarKinematicViscosity(rho, laminarMu);
    const double sqrtK = std::sqrt(std::max(kVal, 0.0));
    const double term1 = 2.0 * sqrtK / std::max(c.betaStar * omega * d, 1.0e-30);
    const double term2 = 500.0 * nu / std::max(d * d * omega, 1.0e-30);
    const double arg2 = std::max(term1, term2);
    return std::tanh(std::pow(std::max(arg2, 0.0), 2.0));
}


} // namespace

void KOmegaSSTModel::initialize(const Field& flow,
                                ScalarFields& state,
                                const FDM::TurbulenceConfig& config) {
    initializeScalar(flow, state, config.scalars.kInitial, ScalarSlot::K);
    initializeScalar(flow, state, config.scalars.omegaInitial, ScalarSlot::Omega);
}

void KOmegaSSTModel::applyPhysicalBoundary(const Field& flow,
                                   ScalarFields& state,
                                   const FDM::TurbulenceConfig& config) {

    int ng = flow.NG();
    int nx = flow.NX();
    int ny = flow.NY();
    int nz = flow.NZ();
    for (int k = ng; k < nz + ng; ++k) {
        for (int j = ng; j < ny + ng; ++j) {
            for (int i = ng; i < nx + ng; ++i) {
                if (wall_ && flow.CellFlag(i,j,k)!=FLUID_CELL) continue;
                state.K(i, j, k) = clampPositive(state.K(i, j, k), config.coefficients.kFloor);
                state.Omega(i, j, k) = clampPositive(state.Omega(i, j, k), config.coefficients.omegaFloor);
            }
        }
    }
    if (wall_) wall_->applySST(flow,config.laminarDynamicViscosity,
        state.values(ScalarSlot::K),state.values(ScalarSlot::Omega),state.values(ScalarSlot::EddyMu));
}

FDM::TransportCellClosure KOmegaSSTModel::onceStepClosure(const Field& flow,const ScalarFields& state,
        const FDM::TurbulenceConfig& config,int i,int j,int k) const {
    const auto& coeff = config.coefficients;
    const SSTConstants c;

    double rho   = std::max(flow(i, j, k, RHO), 1.0e-12);
    double kVal  = clampPositive(state.K(i, j, k),     coeff.kFloor);
    double wVal  = clampPositive(state.Omega(i, j, k), coeff.omegaFloor);

    double sMag = strainRateMagnitude(flow, i, j, k);
    double muT  = eddyDynamicViscosity(flow, state, config, i, j, k);
    double f1 = blendingF1(
        flow, state, i, j, k, rho, kVal, wVal,
        config.laminarDynamicViscosity, c,wall_!=nullptr);
    double sigmaK = f1 * c.sigmaK1 + (1.0 - f1) * c.sigmaK2;
    double sigmaW = f1 * c.sigmaW1 + (1.0 - f1) * c.sigmaW2;
    double beta = f1 * c.beta1 + (1.0 - f1) * c.beta2;
    double gamma = f1 * c.gamma1 + (1.0 - f1) * c.gamma2;

    double production = std::min(muT * sMag * sMag,
                                 20.0 * c.betaStar * rho * kVal * wVal);

    double dk = production / rho
              - c.betaStar * kVal * wVal;
    double domega = gamma * production / std::max(muT, 1.0e-30)
                  - beta * wVal * wVal;

    const Vector3 gradK = scalarPhysicalGradient(
        flow, state, ScalarSlot::K, i, j, k);
    const Vector3 gradW = scalarPhysicalGradient(
        flow, state, ScalarSlot::Omega, i, j, k);
    double crossDiff = 2.0 * (1.0 - f1) * rho * c.sigmaW2
                     * SF::dot(gradK, gradW)
                     / std::max(wVal, coeff.omegaFloor);
    FDM::TransportCellClosure result;result.density=rho;
    result.values={kVal,wVal};result.source={dk,domega};result.diffusionWeight={sigmaK,sigmaW};result.crossDiffusion=crossDiff;
    return result;
}

double KOmegaSSTModel::eddyDynamicViscosity(const Field& flow,
                                            const ScalarFields& state,
                                            const FDM::TurbulenceConfig& config,
                                            int i, int j, int k) const {
    const SSTConstants c;
    double rho = std::max(flow(i, j, k, RHO), 1.0e-12);
    double kVal = clampPositive(state.clamped(ScalarSlot::K, i, j, k), config.coefficients.kFloor);
    double omegaVal = clampPositive(state.clamped(ScalarSlot::Omega, i, j, k), config.coefficients.omegaFloor);
    double f2 = blendingF2(
        flow, i, j, k, rho, kVal, omegaVal,
        config.laminarDynamicViscosity, c,wall_!=nullptr);
    double omegaMag = vorticityMagnitude(flow, i, j, k);
    return rho * c.a1 * kVal
         / std::max(c.a1 * omegaVal, omegaMag * f2);
}

FDM::TransportStageCoefficients KOmegaSSTModel::stageCoefficients(const Field& flow,
        const FDM::TransportStageSamples& samples,const FDM::TurbulenceConfig& config,int i,int j,int k) const {
    const SSTConstants c;
    const double rho=flow(i,j,k,RHO),kv=samples[0](i,j,k),wv=samples[1](i,j,k);
    if(!std::isfinite(rho) || rho<=0 || !std::isfinite(kv) || kv<=0 || !std::isfinite(wv) || wv<=0)
        throw std::runtime_error("Stage SST requires positive finite rho/k/omega; no positivity repair.");
    const auto gradK=CENTRAL2::gradSampled(flow,samples[0],i,j,k);
    const auto gradW=CENTRAL2::gradSampled(flow,samples[1],i,j,k);
    const double product=dot(gradK,gradW),d=wallDistance(flow,i,j,k,wall_!=nullptr);
    const double nu=laminarKinematicViscosity(rho,config.laminarDynamicViscosity);
    const double cd=std::max(2*rho*c.sigmaW2*product/std::max(wv,1e-30),1e-20);
    const double t1=std::sqrt(kv)/std::max(c.betaStar*wv*d,1e-30);
    const double t2=500*nu/std::max(d*d*wv,1e-30);
    const double t3=4*rho*c.sigmaW2*kv/std::max(cd*d*d,1e-30);
    FDM::TransportStageCoefficients value;
    value.blending1=std::tanh(std::pow(std::max(std::min(std::max(t1,t2),t3),0.),4.));
    value.blending2=std::tanh(std::pow(std::max(std::max(2*t1,t2),0.),2.));
    value.dynamicViscosity=rho*c.a1*kv/std::max(c.a1*wv,vorticityMagnitude(flow,i,j,k)*value.blending2);
    const double f1=value.blending1,s=strainRateMagnitude(flow,i,j,k);
    const double production=std::min(value.dynamicViscosity*s*s,20*c.betaStar*rho*kv*wv);
    value.diffusivity={(config.laminarDynamicViscosity+(f1*c.sigmaK1+(1-f1)*c.sigmaK2)*value.dynamicViscosity)/rho,
        (config.laminarDynamicViscosity+(f1*c.sigmaW1+(1-f1)*c.sigmaW2)*value.dynamicViscosity)/rho};
    value.production={production,rho*(f1*c.gamma1+(1-f1)*c.gamma2)*production/std::max(value.dynamicViscosity,1e-30)};
    value.sink={c.betaStar*rho*wv,(f1*c.beta1+(1-f1)*c.beta2)*rho*wv};
    value.crossSource=2*(1-f1)*rho*c.sigmaW2*product/std::max(wv,config.coefficients.omegaFloor);
    return value;
}

} // namespace RAS
} // namespace Turbulence
} // namespace SF
