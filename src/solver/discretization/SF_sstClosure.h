#pragma once
/// @brief SST 窄数学核；不拥有模型、物理数组或时间循环。
#include "methods/numerics/structured/SF_vectorCalculus.h"
#include "methods/numerics/viscous/SF_viscous.h"
#include <functional>
namespace SF::Discretization {
/// @brief 从 WHAT 操作数冻结的有限模型参数。
struct SstParameters {
    double molecularViscosity=0,omegaFloor=0,betaStar=0,sigmaW2=0,gamma1=0,gamma2=0,a1=0;
};
/// @brief 一次共享计算的五个公开结果及局部 F2。
struct SstValues {
    double dynamicViscosity=0,blending1=0,crossSource=0;
    std::array<double,2> production{};
};
using SstSamples=std::array<std::function<double(int,int,int)>,2>;
/// @brief 独立黏度计算只消费本式所需数据，不发布 F1/production/cross。
inline double sstEddyViscosity(double rho,double k,double omega,double distance,
        const Tensor3& gradient,const SstParameters& c) {
    if(!std::isfinite(rho) || rho<=0 || !std::isfinite(k) || k<=0 || !std::isfinite(omega) || omega<=0)
        throw std::runtime_error("SST viscosity needs positive finite density/k/omega");
    const double nu=std::max(c.molecularViscosity,0.)/std::max(rho,1e-12);
    const double t1=std::sqrt(k)/std::max(c.betaStar*omega*distance,1e-30);
    const double t2=500*nu/std::max(distance*distance*omega,1e-30);
    const double f2=std::tanh(std::pow(std::max(std::max(2*t1,t2),0.),2.));
    return rho*c.a1*k/std::max(c.a1*omega,norm(Math::curl(gradient))*f2);
}
/// @brief packed flow 是已声明的 rho/rhoU backing；distance 是显式几何端口。
inline SstValues sstCentral2(const Field& flow,const SstSamples& samples,
        const SstParameters& c,const std::function<double(int,int,int)>& density,
        const std::function<Tensor3(int,int,int)>& velocityGradient,
        const std::function<double(int,int,int)>& distance,double publishedViscosity,int i,int j,int k) {
    const double rho=density(i,j,k),kv=samples[0](i,j,k),wv=samples[1](i,j,k);
    if(!std::isfinite(rho) || rho<=0 || !std::isfinite(kv) || kv<=0 || !std::isfinite(wv) || wv<=0)
        throw std::runtime_error("Stage SST requires positive finite rho/k/omega; no positivity repair.");
    const auto gradK=CENTRAL2::gradSampled(flow,samples[0],i,j,k);
    const auto gradW=CENTRAL2::gradSampled(flow,samples[1],i,j,k);
    const double product=dot(gradK,gradW),d=distance(i,j,k);
    const double nu=std::max(c.molecularViscosity,0.)/std::max(rho,1e-12);
    const double cd=std::max(2*rho*c.sigmaW2*product/std::max(wv,1e-30),1e-20);
    const double t1=std::sqrt(kv)/std::max(c.betaStar*wv*d,1e-30);
    const double t2=500*nu/std::max(d*d*wv,1e-30);
    const double t3=4*rho*c.sigmaW2*kv/std::max(cd*d*d,1e-30);
    SstValues value;
    value.blending1=std::tanh(std::pow(std::max(std::min(std::max(t1,t2),t3),0.),4.));
    value.dynamicViscosity=publishedViscosity;
    const auto strain=Math::symm(velocityGradient(i,j,k));
    const double f1=value.blending1,s=std::sqrt(std::max(2.0*Math::doubleDot(strain,strain),0.0));
    const double production=std::min(value.dynamicViscosity*s*s,20*c.betaStar*rho*kv*wv);
    value.production={production,rho*(f1*c.gamma1+(1-f1)*c.gamma2)*production/std::max(value.dynamicViscosity,1e-30)};
    value.crossSource=2*(1-f1)*rho*c.sigmaW2*product/std::max(wv,c.omegaFloor);
    return value;
}

} // namespace SF::Discretization
