#include "solver/algorithm/pressure/SF_fixedTimeMath.h"
#include "core/config/types/SF_pressureConfigTypes.h"

#include <cmath>
#include <stdexcept>

int main() {
    // (U-U_n)/dt + a*U_k=b. The second outer iterate must retain U_n=2,
    // even after the first candidate has changed the working iterate.
    constexpr double base=2.0,dt=0.5,a=1.0,b=4.0;
    const double first=SF::Pressure::fixedTimeCandidate(
        base,dt*(b-a*base));
    const double second=SF::Pressure::fixedTimeCandidate(
        base,dt*(b-a*first));
    const double incorrectlyAdvanced=SF::Pressure::fixedTimeCandidate(
        first,dt*(b-a*first));
    if (std::abs(first-3.0)>1e-14 || std::abs(second-2.5)>1e-14
        || std::abs(incorrectlyAdvanced-3.5)>1e-14
        || std::abs(second-incorrectlyAdvanced)<0.5)
        throw std::runtime_error("fixed-time temporal base changed across outer iterations");
    if (std::abs(SF::Pressure::relaxSolution(2.0,3.0,0.6)-2.6)>1e-14
        || SF::Pressure::relaxSolution(2.0,3.0,1.0)!=3.0)
        throw std::runtime_error("fixed-point solution relaxation changed");
    // An iterated spatial residual is not forward Euler, even with one clock commit.
    double decay=base;
    for (int k=0;k<100;++k)
        decay=SF::Pressure::fixedTimeCandidate(base,-dt*a*decay);
    if (std::abs(decay-base/(1.0+dt*a))>1e-13
        || std::abs(decay-base*(1.0-dt*a))<0.1)
        throw std::runtime_error("fixed-time residual time-level contract changed");
    SF::FDM::PressureCorrectionConfig config;
    config.reference.referenceCell=0;
    config.reference.referencePressure=0.0;
    SF::FDM::validatePressureCorrectionConfig(config);
    config.coupling.preset=SF::FDM::PressureCouplingPreset::PISO;
    config.coupling.outerCorrectors=2;
    bool invalidPiso=false;
    try { SF::FDM::validatePressureCorrectionConfig(config); }
    catch (const std::invalid_argument&) { invalidPiso=true; }
    if (!invalidPiso) throw std::runtime_error("PISO accepted two outer iterations");
    config.coupling.preset=SF::FDM::PressureCouplingPreset::SIMPLE;
    config.coupling.outerCorrectors=2;
    config.coupling.pressureCorrectors=2;
    bool invalidSimple=false;
    try { SF::FDM::validatePressureCorrectionConfig(config); }
    catch (const std::invalid_argument&) { invalidSimple=true; }
    if (!invalidSimple) throw std::runtime_error("SIMPLE accepted two pressure corrections");
}
