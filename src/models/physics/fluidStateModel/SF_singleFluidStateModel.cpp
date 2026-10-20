/// @file SF_singleFluidStateModel.cpp
/// @brief 守恒变量 FluidStateModel 与工厂实现。

#include "SF_singleFluidStateModel.h"
#include "models/physics/EOS/SF_perfectGasEOS.h"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace SF::Physics::FluidStateModel {

SingleFluidStateModel::SingleFluidStateModel(
        EOS::ModelPtr eos,
        double dynamicViscosity,
        double thermalConductivity)
    : eos_(std::move(eos)),
      dynamicViscosity_(dynamicViscosity),
      thermalConductivity_(thermalConductivity) {
    if (!eos_) {
        throw std::runtime_error(
            "SingleFluidStateModel requires EOS.");
    }
    if (!std::isfinite(dynamicViscosity_)||dynamicViscosity_<0.0
        ||!std::isfinite(thermalConductivity_)
        ||thermalConductivity_<0.0) {
        throw std::runtime_error(
            "SingleFluidStateModel transport properties are invalid.");
    }
    variables_={{"rho"},{"rhoU"},{"rhoV"},{"rhoW"},{"rhoE"}};
    layout_=State::StateLayout({
        {"rho",State::FieldLocation::Cell,
         State::ConservationKind::Conservative,
         State::UpdatePolicy::Explicit,{},"kg/m3"},
        {"rhoU",State::FieldLocation::Cell,
         State::ConservationKind::Conservative,
         State::UpdatePolicy::Explicit,{},"kg/(m2 s)"},
        {"rhoV",State::FieldLocation::Cell,
         State::ConservationKind::Conservative,
         State::UpdatePolicy::Explicit,{},"kg/(m2 s)"},
        {"rhoW",State::FieldLocation::Cell,
         State::ConservationKind::Conservative,
         State::UpdatePolicy::Explicit,{},"kg/(m2 s)"},
        {"rhoE",State::FieldLocation::Cell,
         State::ConservationKind::Conservative,
         State::UpdatePolicy::Explicit,{},"J/m3"}});
}

ThermodynamicState SingleFluidStateModel::close(
        const double* q,int count) const {
    requireEnergyInput(q,count);
    const double rho=q[0];
    ThermodynamicState result;
    result.density=rho;
    for (int direction=0;direction<3;++direction) {
        result.velocity[(size_t)direction]=q[1+direction]/rho;
    }
    const double kinetic=0.5*rho
        *(result.velocity[0]*result.velocity[0]
          +result.velocity[1]*result.velocity[1]
          +result.velocity[2]*result.velocity[2]);
    result.internalEnergyDensity=q[4]-kinetic;
    const EOS::State thermo=eos_->fromDensityEnergy(
        rho,result.internalEnergyDensity/rho);
    result.pressure=thermo.pressure;
    result.temperature=thermo.temperature;
    result.soundSpeed=thermo.soundSpeed;
    result.dynamicViscosity=dynamicViscosity_;
    result.thermalConductivity=thermalConductivity_;
    result.phaseMass={rho};
    result.volumeFraction={1.0};
    return result;
}

double SingleFluidStateModel::totalEnergyFromPressure(
        const double* q,int count,double pressure) const {
    requireEnergyInput(q,count);
    return totalEnergyFromThermo(
        q,count,eos_->fromDensityPressure(q[0],pressure));
}

double SingleFluidStateModel::totalEnergyFromTemperature(
        const double* q,int count,double temperature) const {
    requireEnergyInput(q,count);
    return totalEnergyFromThermo(
        q,count,eos_->fromDensityTemperature(q[0],temperature));
}

std::optional<double> SingleFluidStateModel::perfectGasGamma() const {
    const auto* perfectGas = dynamic_cast<const EOS::PerfectGasEOS*>(eos_.get());
    if (!perfectGas) return std::nullopt;
    return perfectGas->gamma();
}

double SingleFluidStateModel::totalEnergyFromThermo(
        const double* q,int count,const EOS::State& thermo) {
    requireEnergyInput(q,count);
    const double inverseDensity=1.0/q[0];
    const double kinetic=0.5*inverseDensity
        *(q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
    return q[0]*thermo.internalEnergy+kinetic;
}

void SingleFluidStateModel::requireEnergyInput(
        const double* q,int count) {
    if (!q||count!=5||!std::isfinite(q[0])||q[0]<=0.0) {
        throw std::runtime_error(
            "SingleFluidStateModel energy inversion received invalid Q.");
    }
    for (int component=1;component<4;++component) {
        if (!std::isfinite(q[component])) {
            throw std::runtime_error(
                "SingleFluidStateModel received non-finite momentum.");
        }
    }
}

} // namespace SF::Physics::FluidStateModel
