/// @file SF_interphase.cpp
/// @brief 双欧拉相间力、热或质量传递闭式模型实现。

#include "SF_interphase.h"

#include <cmath>
#include <sstream>
#include <stdexcept>

namespace SF::Physics::PhaseSystems::Interphase {
namespace {

size_t phaseIndex(const Context& context, const std::string& name) {
    for (size_t phase = 0; phase < context.phases.size(); ++phase) {
        if (Multiphase::normalizePhaseName(context.phases[phase].name)
            == Multiphase::normalizePhaseName(name)) {
            return phase;
        }
    }
    throw std::runtime_error(
        "Interphase phase pair references unknown phase '" + name + "'.");
}

} // namespace

void addConservativePairForce(
        const PairContext& context,
        PhaseEquationSources& sources,
        int cell,
        const std::array<double, 3>& forceOnContinuous) {
    const auto& continuous = context.system.phases[context.continuous];
    const auto& dispersed = context.system.phases[context.dispersed];
    double continuousHeating = 0.0;
    double dispersedHeating = 0.0;
    for (int component = 0; component < 3; ++component) {
        const double force = forceOnContinuous[(size_t)component];
        if (!std::isfinite(force)) {
            throw std::runtime_error(
                "Interphase pair force contains a non-finite component.");
        }
        const double interfaceVelocity = 0.5 * (
            continuous.primitive.velocity[(size_t)component]
                .values()[(size_t)cell]
            + dispersed.primitive.velocity[(size_t)component]
                .values()[(size_t)cell]);
        continuousHeating += force * (
            interfaceVelocity
            - continuous.primitive.velocity[(size_t)component]
                .values()[(size_t)cell]);
        dispersedHeating -= force * (
            interfaceVelocity
            - dispersed.primitive.velocity[(size_t)component]
                .values()[(size_t)cell]);
    }
    std::array<double,3> transfer{};
    for(int component=0;component<3;++component) {
        transfer[(size_t)component]=-forceOnContinuous[(size_t)component];
    }
    sources.transferLedger.addInternalMomentumTransfer(
        cell,(int)context.continuous,(int)context.dispersed,transfer);
    sources.transferLedger.addInternalMechanicalEnergy(
        cell,(int)context.continuous,continuousHeating);
    sources.transferLedger.addInternalMechanicalEnergy(
        cell,(int)context.dispersed,dispersedHeating);
    sources.interphaseMechanicalHeating.values()[(size_t)cell]
        += continuousHeating + dispersedHeating;
}

void compute(const Context& context, PhaseEquationSources& sources) {
    if (!std::isfinite(context.dt) || context.dt <= 0.0)
        throw std::runtime_error("Interphase models require finite positive dt.");
    sources.clear();
    if (sources.momentumCouplings.size()
        != context.config.eulerianEulerian.phasePairs.size()) {
        throw std::runtime_error(
            "Interphase phase-pair registry and source workspace differ.");
    }
    for (size_t pairIndex = 0;
         pairIndex < context.config.eulerianEulerian.phasePairs.size();
         ++pairIndex) {
        const auto& options =
            context.config.eulerianEulerian.phasePairs[pairIndex];
        PairContext pair{
            context, options,
            phaseIndex(context, options.continuousPhase),
            phaseIndex(context, options.dispersedPhase),
            pairIndex};
        addDragAndVirtualMass(pair, sources);
        addLift(pair, sources);
        addWallLubrication(pair, sources);
        addTurbulentDispersion(pair, sources);
        addHeatTransfer(pair, sources);
    }
    addPhaseChange(context, sources);
    sources.transferLedger.validate();

    for (int n=0;n<context.geometry.TotalSize();++n) {
        for(size_t phase=0;phase<context.phases.size();++phase) {
            sources.mass[phase].values()[(size_t)n]+=
                sources.transferLedger.phaseMassSource(n,(int)phase);
            sources.energy[phase].values()[(size_t)n]+=
                sources.transferLedger.phaseEnergySource(n,(int)phase);
            for(int component=0;component<3;++component) {
                sources.momentum[phase][(size_t)component]
                    .values()[(size_t)n]+=
                    sources.transferLedger.phaseMomentumSource(
                        n,(int)phase,component);
            }
        }
        double mass=0.0,energy=0.0;
        double massScale=1.0, energyScale=1.0;
        std::array<double,3> momentum{0.0,0.0,0.0};
        std::array<double,3> momentumScale{1.0,1.0,1.0};
        for(size_t k=0;k<context.phases.size();++k){
            const double phaseMass = sources.mass[k].values()[(size_t)n];
            const double phaseEnergy = sources.energy[k].values()[(size_t)n];
            mass += phaseMass;
            energy += phaseEnergy;
            massScale += std::abs(phaseMass);
            energyScale += std::abs(phaseEnergy);
            for(int d=0;d<3;++d) {
                const double phaseMomentum =
                    sources.momentum[k][(size_t)d].values()[(size_t)n];
                momentum[(size_t)d] += phaseMomentum;
                momentumScale[(size_t)d] += std::abs(phaseMomentum);
            }
        }
        const double dissipation =
            sources.interphaseMechanicalHeating.values()[(size_t)n];
        energyScale += std::abs(dissipation);
        const bool massFailed = std::abs(mass) > 1.0e-12 * massScale;
        const bool energyFailed =
            std::abs(energy - dissipation) > 1.0e-12 * energyScale;
        bool momentumFailed = false;
        for (int d = 0; d < 3; ++d) {
            momentumFailed = momentumFailed
                || std::abs(momentum[(size_t)d])
                    > 1.0e-12 * momentumScale[(size_t)d];
        }
        if(massFailed || energyFailed || momentumFailed){
            std::ostringstream message;
            message.precision(17);
            message
                << "Interphase source ledger is not pairwise conservative"
                << " at cell " << n
                << ": massImbalance=" << mass
                << ", massScale=" << massScale
                << ", momentumImbalance=(" << momentum[0] << ", "
                << momentum[1] << ", " << momentum[2] << ")"
                << ", momentumScale=(" << momentumScale[0] << ", "
                << momentumScale[1] << ", " << momentumScale[2] << ")"
                << ", energyImbalance=" << energy - dissipation
                << ", energyScale=" << energyScale << '.';
            throw std::runtime_error(message.str());
        }
    }
}

} // namespace SF::Physics::PhaseSystems::Interphase
