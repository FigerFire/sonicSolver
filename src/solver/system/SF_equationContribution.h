#pragma once

/// @file SF_equationContribution.h
/// @brief 模块向现有 ResolvedSimulationSystem 声明数学贡献的窄接口。

#include "SF_resolvedSimulationSystem.h"

namespace SF::System {

class EquationSystemBuilder {
public:
    explicit EquationSystemBuilder(ResolvedSimulationSystem& system)
        : system_(system) {}

    void addUnknown(UnknownDescriptor unknown);
    void addEquation(
        EquationDescriptor descriptor, Equation::Definition definition);
    void addTerm(const std::string& equationId, Equation::Term term);
    void addConstraint(ConstraintDescriptor constraint);
    void addSolveBlock(SolveBlock block);
    void addClosure(std::string closure);
    void require(ExecutionRequirement requirement);

private:
    ResolvedSimulationSystem& system_;
};

/// @brief Contribution 只描述数学增量，不拥有或调用 timestep lifecycle。
class IEquationContribution {
public:
    virtual ~IEquationContribution() = default;
    virtual void contribute(EquationSystemBuilder& system) const = 0;
};

} // namespace SF::System
