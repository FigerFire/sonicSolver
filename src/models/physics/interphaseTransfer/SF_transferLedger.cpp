/// @file SF_transferLedger.cpp
/// @brief 相间质量、动量与能量守恒账本实现。

#include "SF_transferLedger.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace SF::Physics::InterphaseTransfer {

SourceBundle::SourceBundle(int cellCount, int variableCount)
    : cellCount_(cellCount), variableCount_(variableCount) {
    if (cellCount_ <= 0 || variableCount_ <= 0) {
        throw std::runtime_error("SourceBundle dimensions must be positive.");
    }
    values_.assign((size_t)cellCount_*(size_t)variableCount_, 0.0);
}

double& SourceBundle::operator()(int cell, int variable) {
    if (cell < 0 || cell >= cellCount_ || variable < 0
        || variable >= variableCount_) {
        throw std::runtime_error("SourceBundle index is out of range.");
    }
    return values_[(size_t)cell*(size_t)variableCount_+(size_t)variable];
}

double SourceBundle::operator()(int cell, int variable) const {
    if (cell < 0 || cell >= cellCount_ || variable < 0
        || variable >= variableCount_) {
        throw std::runtime_error("SourceBundle index is out of range.");
    }
    return values_[(size_t)cell*(size_t)variableCount_+(size_t)variable];
}

TransferLedger::TransferLedger(int cellCount, int phaseCount)
    : cellCount_(cellCount), phaseCount_(phaseCount) {
    if (cellCount_ <= 0 || phaseCount_ <= 0) {
        throw std::runtime_error("TransferLedger dimensions must be positive.");
    }
    const size_t phaseValues = (size_t)cellCount_*(size_t)phaseCount_;
    mass_.assign(phaseValues, 0.0);
    momentum_.assign(phaseValues*3, 0.0);
    energy_.assign(phaseValues, 0.0);
    mixtureEnergy_.assign((size_t)cellCount_, 0.0);
    wallEnergy_.assign((size_t)cellCount_, 0.0);
    externalEnergy_.assign((size_t)cellCount_, 0.0);
    mechanicalEnergy_.assign((size_t)cellCount_,0.0);
    phaseExternalEnergy_.assign((size_t)cellCount_,0.0);
}

void TransferLedger::requireCell(int cell) const {
    if (cell < 0 || cell >= cellCount_) {
        throw std::runtime_error("TransferLedger cell index is out of range.");
    }
}

void TransferLedger::requirePhase(int phase) const {
    if (phase < 0 || phase >= phaseCount_) {
        throw std::runtime_error("TransferLedger phase index is out of range.");
    }
}

size_t TransferLedger::phaseIndex(int cell, int phase) const {
    requireCell(cell);
    requirePhase(phase);
    return (size_t)cell*(size_t)phaseCount_+(size_t)phase;
}

size_t TransferLedger::momentumIndex(
        int cell, int phase, int component) const {
    if (component < 0 || component >= 3) {
        throw std::runtime_error(
            "TransferLedger momentum component is out of range.");
    }
    return phaseIndex(cell, phase)*3+(size_t)component;
}

void TransferLedger::addInternalMassTransfer(
        int cell, int donorPhase, int receiverPhase, double rate) {
    if (donorPhase == receiverPhase) {
        throw std::runtime_error(
            "TransferLedger internal mass transfer needs distinct phases.");
    }
    if (!std::isfinite(rate) || rate < 0.0) {
        throw std::runtime_error(
            "TransferLedger internal mass-transfer rate must be finite and non-negative.");
    }
    mass_[phaseIndex(cell, donorPhase)] -= rate;
    mass_[phaseIndex(cell, receiverPhase)] += rate;
}

void TransferLedger::addInternalMomentumTransfer(
        int cell, int donorPhase, int receiverPhase,
        const std::array<double, 3>& rate) {
    if (donorPhase == receiverPhase) {
        throw std::runtime_error(
            "TransferLedger internal momentum transfer needs distinct phases.");
    }
    for (int d = 0; d < 3; ++d) {
        if (!std::isfinite(rate[(size_t)d])) {
            throw std::runtime_error(
                "TransferLedger internal momentum rate is non-finite.");
        }
        momentum_[momentumIndex(cell, donorPhase, d)] -= rate[(size_t)d];
        momentum_[momentumIndex(cell, receiverPhase, d)] += rate[(size_t)d];
    }
}

void TransferLedger::addInternalEnergyTransfer(
        int cell, int donorPhase, int receiverPhase, double rate) {
    if (donorPhase == receiverPhase) {
        throw std::runtime_error(
            "TransferLedger internal energy transfer needs distinct phases.");
    }
    if (!std::isfinite(rate)) {
        throw std::runtime_error(
            "TransferLedger internal energy rate is non-finite.");
    }
    energy_[phaseIndex(cell, donorPhase)] -= rate;
    energy_[phaseIndex(cell, receiverPhase)] += rate;
}

void TransferLedger::addInternalMechanicalEnergy(
        int cell,int phase,double rate) {
    if(!std::isfinite(rate)) {
        throw std::runtime_error(
            "TransferLedger mechanical energy conversion must be finite.");
    }
    energy_[phaseIndex(cell,phase)]+=rate;
    mechanicalEnergy_[(size_t)cell]+=rate;
}

void TransferLedger::addWallEnergyInput(int cell, double rate) {
    requireCell(cell);
    if (!std::isfinite(rate)) {
        throw std::runtime_error("TransferLedger wall-energy input is non-finite.");
    }
    wallEnergy_[(size_t)cell] += rate;
    mixtureEnergy_[(size_t)cell] += rate;
}

void TransferLedger::addWallEnergyInput(
        int cell,int phase,double rate) {
    addWallEnergyInput(cell,rate);
    energy_[phaseIndex(cell,phase)]+=rate;
    phaseExternalEnergy_[(size_t)cell]+=rate;
}

void TransferLedger::addExternalEnergyInput(int cell, double rate) {
    requireCell(cell);
    if (!std::isfinite(rate)) {
        throw std::runtime_error("TransferLedger external-energy input is non-finite.");
    }
    externalEnergy_[(size_t)cell] += rate;
    mixtureEnergy_[(size_t)cell] += rate;
}

double TransferLedger::phaseMassSource(int cell, int phase) const {
    return mass_[phaseIndex(cell, phase)];
}

double TransferLedger::phaseMomentumSource(
        int cell, int phase, int component) const {
    return momentum_[momentumIndex(cell, phase, component)];
}

double TransferLedger::phaseEnergySource(int cell, int phase) const {
    return energy_[phaseIndex(cell, phase)];
}

double TransferLedger::mixtureEnergySource(int cell) const {
    requireCell(cell);
    return mixtureEnergy_[(size_t)cell];
}

double TransferLedger::wallEnergyInput(int cell) const {
    requireCell(cell);
    return wallEnergy_[(size_t)cell];
}

double TransferLedger::externalEnergyInput(int cell) const {
    requireCell(cell);
    return externalEnergy_[(size_t)cell];
}

double TransferLedger::mechanicalEnergyConversion(int cell) const {
    requireCell(cell);
    return mechanicalEnergy_[(size_t)cell];
}

void TransferLedger::validate(double relativeTolerance) const {
    if (!std::isfinite(relativeTolerance) || relativeTolerance < 0.0) {
        throw std::runtime_error(
            "TransferLedger validation tolerance must be finite and non-negative.");
    }
    for (int cell = 0; cell < cellCount_; ++cell) {
        double massSum = 0.0;
        double massScale = 0.0;
        std::array<double, 3> momentumSum{0.0, 0.0, 0.0};
        std::array<double, 3> momentumScale{0.0, 0.0, 0.0};
        double energySum = 0.0;
        double energyScale = 0.0;
        for (int phase = 0; phase < phaseCount_; ++phase) {
            const double mass = phaseMassSource(cell, phase);
            const double energy = phaseEnergySource(cell, phase);
            if (!std::isfinite(mass) || !std::isfinite(energy)) {
                throw std::runtime_error(
                    "TransferLedger contains a non-finite phase source.");
            }
            massSum += mass;
            massScale += std::abs(mass);
            energySum += energy;
            energyScale += std::abs(energy);
            for (int d = 0; d < 3; ++d) {
                const double momentum = phaseMomentumSource(cell, phase, d);
                if (!std::isfinite(momentum)) {
                    throw std::runtime_error(
                        "TransferLedger contains a non-finite momentum source.");
                }
                momentumSum[(size_t)d] += momentum;
                momentumScale[(size_t)d] += std::abs(momentum);
            }
        }
        const double expectedEnergy=
            mechanicalEnergyConversion(cell)
            +phaseExternalEnergy_[(size_t)cell];
        energyScale+=std::abs(expectedEnergy);
        if (std::abs(massSum) > relativeTolerance*std::max(massScale, 1.0)
            || std::abs(energySum-expectedEnergy)
                >relativeTolerance*std::max(energyScale,1.0)) {
            throw std::runtime_error(
                "TransferLedger internal mass/energy exchange is not conservative at cell "
                + std::to_string(cell) + ".");
        }
        for (int d = 0; d < 3; ++d) {
            if (std::abs(momentumSum[(size_t)d])
                > relativeTolerance*std::max(momentumScale[(size_t)d], 1.0)) {
                throw std::runtime_error(
                    "TransferLedger internal momentum exchange is not conservative at cell "
                    + std::to_string(cell) + ".");
            }
        }
        const double external = wallEnergyInput(cell)+externalEnergyInput(cell);
        const double scale = std::max({std::abs(external),
                                       std::abs(mixtureEnergySource(cell)), 1.0});
        if (std::abs(mixtureEnergySource(cell)-external)
            > relativeTolerance*scale) {
            throw std::runtime_error(
                "TransferLedger external mixture-energy budget is not closed at cell "
                + std::to_string(cell) + ".");
        }
    }
}

} // namespace SF::Physics::InterphaseTransfer
