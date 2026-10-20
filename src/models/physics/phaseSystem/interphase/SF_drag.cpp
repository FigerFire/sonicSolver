/// @file SF_drag.cpp
/// @brief 双欧拉相间力、热或质量传递闭式模型实现。

#include "SF_interphase.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace SF::Physics::PhaseSystems::Interphase {

void addDragAndVirtualMass(
        const PairContext& pair,
        PhaseEquationSources& sources) {
    const Context& c = pair.system;
    const auto& opt = pair.options;
    if(Multiphase::normalizeModelType(opt.dragModel)!="schillernaumann")
        throw std::runtime_error("Eulerian drag supports only SchillerNaumann.");
    if(Multiphase::normalizeModelType(opt.virtualMassModel)!="constantcoefficient")
        throw std::runtime_error("Eulerian virtualMass supports only constantCoefficient.");
    const auto& pc=c.properties[pair.continuous];
    const auto& continuous=c.phases[pair.continuous];
    const auto& dispersed=c.phases[pair.dispersed];
    const double configuredDiameter =
        Multiphase::phasePairDiameter(c.config, opt);
    auto& coupling = sources.momentumCouplings[pair.couplingIndex];
    for(int n=0;n<c.geometry.TotalSize();++n){
        const double ac=continuous.primitive.alpha.values()[(size_t)n];
        const double ad=dispersed.primitive.alpha.values()[(size_t)n];
        const double rho=continuous.primitive.density.values()[(size_t)n];
        const double wallDiameter =
            sources.wallBoilingDepartureDiameter.values()[(size_t)n];
        const double diameter = wallDiameter > 0.0
            ? wallDiameter : configuredDiameter;
        if (!std::isfinite(diameter) || diameter <= 0.0)
            throw std::runtime_error(
                "Eulerian drag requires a positive local bubble diameter.");
        double ur2=0.0;
        for(int d=0;d<3;++d){
            const double ur=dispersed.primitive.velocity[(size_t)d].values()[(size_t)n]
                           -continuous.primitive.velocity[(size_t)d].values()[(size_t)n];
            ur2+=ur*ur;
        }
        const double urMag=std::sqrt(ur2);
        const double temperature = continuous.primitive.temperature
            .values()[(size_t)n];
        const double viscosity =
            Multiphase::phaseViscosity(pc, temperature);
        const double re=rho*urMag*diameter/viscosity;
        double beta=0.0;
        if(re==0.0){
            // Schiller–Naumann 在 Re->0 的解析 Stokes 极限，不是数值兜底。
            beta=18.0*viscosity*ac*ad/
                 (diameter*diameter);
        }else{
            const double cd=re<1000.0
                ? 24.0/re*(1.0+0.15*std::pow(re,0.687)) : 0.44;
            beta=0.75*cd*rho*ac*ad*urMag/diameter;
        }
        const double vm=opt.virtualMassCoefficient*rho*ad/c.dt;
        const double coefficient=beta+vm;
        if(!std::isfinite(coefficient)||coefficient<0.0)
            throw std::runtime_error("Invalid drag/virtual-mass coefficient.");
        coupling.coefficient.values()[(size_t)n]=coefficient;
        coupling.dragCoefficient.values()[(size_t)n]=beta;
        std::array<double, 3> force{};
        for(int d=0;d<3;++d){
            const double uc=continuous.primitive.velocity[(size_t)d].values()[(size_t)n];
            const double ud=dispersed.primitive.velocity[(size_t)d].values()[(size_t)n];
            const double prevUc=c.previousVelocity[pair.continuous][(size_t)d]
                                    .values()[(size_t)n];
            const double prevUd=c.previousVelocity[pair.dispersed][(size_t)d]
                                    .values()[(size_t)n];
            const double drag=beta*(ud-uc);
            const double virtualMass=opt.virtualMassCoefficient*rho*ad
                *((ud-prevUd)-(uc-prevUc))/c.dt;
            force[(size_t)d]=drag+virtualMass;
        }
        addConservativePairForce(pair, sources, n, force);
    }
}

} // namespace SF::Physics::PhaseSystems::Interphase
