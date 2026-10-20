#pragma once

/// @file SF_interphase.h
/// @brief Eulerian–Eulerian 相间力、换热和相变调度。

#include "SF_phaseSystem.h"

namespace SF::Physics::PhaseSystems::Interphase {

/// @brief 所有闭式模型共享的只读输入。
struct Context {
    const Field& geometry;
    const Multiphase::MultiPhaseConfig& config;
    const std::vector<Multiphase::PhaseProperties>& properties;
    const std::vector<PhaseState>& phases;
    const std::vector<PhaseVectorField>& previousVelocity;
    double dt = 0.0;
};

/// @brief 单个 phase-pair 闭式的解析后上下文。
struct PairContext {
    const Context& system;
    const Multiphase::PhasePairModelOptions& options;
    size_t continuous = 0;
    size_t dispersed = 1;
    size_t couplingIndex = 0;
};

/// @brief 计算 drag/lift/virtual-mass/heat/mass-transfer 源。
void compute(const Context& context, PhaseEquationSources& sources);

void addDragAndVirtualMass(const PairContext& context,
                           PhaseEquationSources& sources);
void addLift(const PairContext& context, PhaseEquationSources& sources);
void addWallLubrication(const PairContext& context,
                        PhaseEquationSources& sources);
void addTurbulentDispersion(const PairContext& context,
                            PhaseEquationSources& sources);
void addHeatTransfer(const PairContext& context,
                     PhaseEquationSources& sources);
void addPhaseChange(const Context& context, PhaseEquationSources& sources);

/// @brief 累加一个作用在连续相上的成对力及其公共界面速度功。
void addConservativePairForce(
    const PairContext& context,
    PhaseEquationSources& sources,
    int cell,
    const std::array<double, 3>& forceOnContinuous);

} // namespace SF::Physics::PhaseSystems::Interphase
