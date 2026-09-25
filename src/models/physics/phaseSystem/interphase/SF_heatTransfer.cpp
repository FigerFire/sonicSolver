/// @file SF_heatTransfer.cpp
/// @brief 双欧拉相间力、热或质量传递闭式模型实现。

#include "SF_interphase.h"

#include <cmath>
#include <stdexcept>

namespace SF::Physics::PhaseSystems::Interphase {

void addHeatTransfer(
        const PairContext& pair,
        PhaseEquationSources& sources) {
    const Context& c = pair.system;
    const auto& opt = pair.options;
    if(Multiphase::normalizeModelType(opt.heatTransferModel)!="ranzmarshall")
        throw std::runtime_error("Eulerian heat transfer supports only RanzMarshall.");
    const auto& prop=c.properties[pair.continuous];
    const auto& pc=c.phases[pair.continuous];
    const auto& pd=c.phases[pair.dispersed];
    const double configuredDiameter =
        Multiphase::phasePairDiameter(c.config, opt);
    for(int n=0;n<c.geometry.TotalSize();++n){
        const double wallDiameter =
            sources.wallBoilingDepartureDiameter.values()[(size_t)n];
        const double diameter = wallDiameter > 0.0
            ? wallDiameter : configuredDiameter;
        if(!std::isfinite(diameter)||diameter<=0.0)
            throw std::runtime_error(
                "RanzMarshall requires positive local bubble diameter.");
        const double continuousTemperature =
            pc.primitive.temperature.values()[(size_t)n];
        const double cpValue =
            Multiphase::phaseSpecificHeat(prop, continuousTemperature);
        const double viscosity =
            Multiphase::phaseViscosity(prop, continuousTemperature);
        const double conductivity =
            Multiphase::phaseThermalConductivity(
                prop, continuousTemperature);
        const double pr=cpValue*viscosity/conductivity;
        if(!std::isfinite(pr)||pr<=0.0)
            throw std::runtime_error(
                "RanzMarshall requires positive Prandtl number.");
        double ur2=0.0;
        for(int d=0;d<3;++d){
            const double ur=pd.primitive.velocity[(size_t)d].values()[(size_t)n]
                           -pc.primitive.velocity[(size_t)d].values()[(size_t)n];
            ur2+=ur*ur;
        }
        const double re=pc.primitive.density.values()[(size_t)n]*std::sqrt(ur2)
                       *diameter/viscosity;
        const double nu=2.0+0.6*std::sqrt(re)*std::cbrt(pr);
        const double areaDensity=6.0*pd.primitive.alpha.values()[(size_t)n]
                                 /diameter;
        const double hA=nu*conductivity/diameter*areaDensity;
        const double q=hA*(pd.primitive.temperature.values()[(size_t)n]
                          -pc.primitive.temperature.values()[(size_t)n]);
        if(!std::isfinite(q))throw std::runtime_error("RanzMarshall heat source is non-finite.");
        sources.transferLedger.addInternalEnergyTransfer(
            n,(int)pair.continuous,(int)pair.dispersed,-q);
    }
}

} // namespace SF::Physics::PhaseSystems::Interphase
