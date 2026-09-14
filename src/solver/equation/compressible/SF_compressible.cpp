/// @file SF_compressible.cpp
/// @brief 单流体可压缩守恒方程项及源项表达式的组装。

#include "solver/equation/compressible/SF_compressible.h"

#include "SF_config.h"
#include "solver/discretization/SF_discretization.h"
#include "SF_field.h"
#include "SF_rusanovEOS.h"

namespace SF::Equation::Compressible {

System::System() {
    const Symbol rho{"rho"};
    const Symbol rhoU{"rhoU"};
    const Symbol rhoE{"rhoE"};
    const Symbol massFlux{"massFlux"};
    const Symbol momentumFlux{"momentumFlux"};
    const Symbol energyFlux{"energyFlux"};
    const Symbol mu{"mu"};
    const Symbol conductivity{"conductivity"};
    const Symbol velocity{"U"};
    const Symbol temperature{"T"};
    const Symbol zero{"zero"};
    const Symbol bodyForces{"bodyForces"};
    const Symbol energySources{"energySources"};

    definition_.add(named(
        "mass", ddt(rho) + div(massFlux) == zero));
    definition_.add(named(
        "momentum",
        ddt(rhoU) + div(momentumFlux) + diffusion(mu, velocity)
            == bodyForces));
    definition_.add(named(
        "energy",
        ddt(rhoE) + div(energyFlux)
            + diffusion(conductivity, temperature)
            == energySources));
}

void System::begin(FluxField& fluxField, Residual& residual) const {
    fluxField.clear();
    residual.clear();
}

void System::convection(
        Field& field, FluxField& fluxField, Residual& residual,
        const AssemblyContext& context) const {
    if (!definition_.contains(TermKind::Divergence)) return;
    if (context.thermodynamics
        && context.convectionThermodynamics
            == ConvectionThermodynamicContract::EquationSetRusanov) {
        Numerics::RusanovEOS::div(field, fluxField, residual,
                                  *context.thermodynamics);
        return;
    }
    double gamma = context.config.numerics.idealGasGamma;
    if (context.thermodynamics) {
        const auto equationGamma = context.thermodynamics->perfectGasGamma();
        if (!equationGamma) {
            throw std::runtime_error(
                "Configured high-order convection requires a PerfectGas EquationSet; "
                "generic EquationSet thermodynamics are not supported by the "
                "current characteristic flux implementation.");
        }
        gamma = *equationGamma;
    }
    divDispatch(field, fluxField, residual,
                context.config.numerics.convection,
                context.config.numerics.flux,
                context.timeStep,
                context.config.numerics.ibmBoundary,
                context.config.numerics.ilwOrder,
                gamma);
}

void System::diffusionAndSources(
        Field& field, Residual& residual, const AssemblyContext& context) const {
    if (definition_.contains(TermKind::Diffusion)) {
        laplacianDispatch(
            field, residual,
            context.config.numerics.viscous,
            context.config.numerics.viscousEnabled || context.transport != nullptr,
            context.config.numerics.dynamicViscosity,
            context.config.numerics.prandtl,
            context.config.numerics.idealGasGamma,
            context.config.numerics.idealGasConstant,
            context.transport);
    }
    if (definition_.contains(TermKind::ExplicitSource)
        || definition_.contains(TermKind::ImplicitSource)) {
        SourceTerm::Sp(field, residual, context.config.sources);
    }
}

void System::assemble(
        Field& field, FluxField& fluxField, Residual& residual,
        const AssemblyContext& context) const {
    begin(fluxField, residual);
    convection(field, fluxField, residual, context);
    diffusionAndSources(field, residual, context);
}

} // namespace SF::Equation::Compressible
