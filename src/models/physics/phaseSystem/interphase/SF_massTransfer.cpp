/// @file SF_massTransfer.cpp
/// @brief 双欧拉相间力、热或质量传递闭式模型实现。

#include "SF_interphase.h"

#include "SF_phaseChange.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace SF::Physics::PhaseSystems::Interphase {

void addPhaseChange(const Context& c, PhaseEquationSources& s){
    if(!c.config.phaseChange.enabled
       || PhaseChange::normalizeModel(c.config.phaseChange.model)=="none")return;
    if(PhaseChange::normalizeModel(c.config.phaseChange.model)=="rpi"){
        // RPI 的质量、动量、显热和潜热由 WallHeat provider 的统一 ledger
        // 一次装配，避免壁面质量源与热流分配分别计算。
        return;
    }
    size_t liquid=c.phases.size(),vapor=c.phases.size();
    for(size_t k=0;k<c.properties.size();++k){
        if(c.properties[k].role==Multiphase::PhaseRole::Liquid)liquid=k;
        if(c.properties[k].role==Multiphase::PhaseRole::Gas)vapor=k;
    }
    if(liquid==c.phases.size()||vapor==c.phases.size())
        throw std::runtime_error("Eulerian phaseChange requires one liquid and one gas phase role.");
    const auto rates=PhaseChange::computeRates(
                                               c.geometry,
                                               c.phases[liquid].primitive.alpha,
                                               c.phases[liquid].primitive.temperature,
                                               c.config,c.dt);
    for(int n=0;n<c.geometry.TotalSize();++n){
        const double candidate=rates.mdot[(size_t)n];
        const size_t donor=candidate>=0.0?liquid:vapor;
        const double available=c.phases[donor].primary.phaseMass
                                   .values()[(size_t)n]/c.dt;
        if(!std::isfinite(available)||available<0.0)
            throw std::runtime_error("Eulerian phaseChange found invalid donor phase mass.");
        if(std::abs(candidate)>available*(1.0+1.0e-12)){
            throw std::runtime_error(
                "Eulerian phaseChange candidate would exhaust donor phaseMass; "
                "reduce maxDeltaT or phaseSourceCFL. candidate="
                + std::to_string(candidate) + ", availableRate="
                + std::to_string(available) + ", cell="
                + std::to_string(n) + ".");
        }
        const double mdot=candidate;
        const size_t source=mdot>=0.0?liquid:vapor;
        const size_t receiver=mdot>=0.0?vapor:liquid;
        const double transferMass=std::abs(mdot);
        s.transferLedger.addInternalMassTransfer(
            n,(int)source,(int)receiver,transferMass);
        std::array<double,3> momentum{};
        for(int d=0;d<3;++d){
            momentum[(size_t)d]=transferMass
                *c.phases[source].primitive.velocity[(size_t)d]
                    .values()[(size_t)n];
        }
        s.transferLedger.addInternalMomentumTransfer(
            n,(int)source,(int)receiver,momentum);
        const double h=c.phases[source].primitive.enthalpy
            .values()[(size_t)n];
        const double transferEnergy=transferMass
            *(h+(mdot>=0.0?c.config.phaseChange.latentHeat:0.0));
        s.transferLedger.addInternalEnergyTransfer(
            n,(int)source,(int)receiver,transferEnergy);
    }
}

} // namespace SF::Physics::PhaseSystems::Interphase
