/// @file SF_moduleGraph.cpp
/// @brief 校验物理模块依赖并生成确定的模块执行图。

#include "SF_moduleGraph.h"

#include "solver/algorithm/immersed/SF_immersedStrategy.h"

#include <sstream>
#include <stdexcept>
#include <utility>

namespace SF::Workflow {
namespace {

using FDM::Capability;
using FDM::CapabilitySet;

ModuleDescriptor equationModule(Kind kind) {
    ModuleDescriptor module;
    if (kind == Kind::Homogeneous) {
        module.name = "HomogeneousEquationSet";
        module.provides = {
            Capability::CanonicalConservativeState,
            Capability::PrimitiveState};
    } else if (kind == Kind::EulerianEulerian) {
        module.name = "EulerianEulerianPhaseSystem";
        module.provides = {Capability::PhaseCanonicalState};
    } else {
        module.name = kind == Kind::OneFluidInterface
            ? "OneFluidEquationSet" : "SingleFluidEquationSet";
        module.provides = {
            Capability::CanonicalConservativeState,
            Capability::PrimitiveState,
            Capability::FiveVariableEulerState,
            Capability::CharacteristicEigenSystem};
    }
    return module;
}

} // namespace

void ModuleGraph::add(ModuleDescriptor module) {
    if (module.name.empty()) {
        throw std::runtime_error("ModuleGraph received an unnamed module.");
    }
    modules_.push_back(std::move(module));
}

FDM::CapabilitySet ModuleGraph::providedCapabilities() const {
    CapabilitySet result;
    for (const auto& module : modules_) result.merge(module.provides);
    return result;
}

void ModuleGraph::validate() const {
    const CapabilitySet available = providedCapabilities();
    for (const auto& module : modules_) {
        for (Capability capability : module.requiresAll.values()) {
            if (!available.contains(capability)) {
                throw std::runtime_error(
                    "module '" + module.name + "' requires missing capability '"
                    + FDM::toString(capability) + "'. Available: ["
                    + available.describe() + "].");
            }
        }
        for (const CapabilitySet& alternatives : module.requiresOneOf) {
            bool satisfied = false;
            for (Capability capability : alternatives.values()) {
                satisfied = satisfied || available.contains(capability);
            }
            if (!satisfied) {
                throw std::runtime_error(
                    "module '" + module.name
                    + "' requires one of [" + alternatives.describe()
                    + "]. Available: [" + available.describe() + "].");
            }
        }
    }
}

std::string ModuleGraph::describe() const {
    std::ostringstream os;
    for (size_t i = 0; i < modules_.size(); ++i) {
        if (i != 0) os << " -> ";
        os << modules_[i].name;
    }
    return os.str();
}

ModuleGraph makeModuleGraph(
        const Request& request,
        const Plan& plan,
        const FDM::SolverConfig& config) {
    ModuleGraph graph;
    graph.add(equationModule(plan.kind));

    ModuleDescriptor flow;
    if (config.numerics.solver == FDM::SolverAlgorithm::PressureBased) {
        flow.name = "PressureBasedAlgorithm";
        flow.provides = {
            Capability::PredictedConservativeState,
            Capability::PressureCorrection,
            Capability::PressureJumpConsumer};
    } else {
        flow.name = "DensityBasedAlgorithm";
        flow.provides.add(Capability::PredictedConservativeState);
    }
    flow.provides.add(Capability::SingleFieldExecution);
    flow.provides.add(Capability::MultiFieldExecution);
    graph.add(std::move(flow));

    ModuleDescriptor geometry;
    geometry.name = request.ibm
        ? "ImmersedBoundaryGeometry" : "BodyFittedGeometry";
    geometry.provides.add(request.ibm
        ? Capability::ImmersedBoundaryGeometry
        : Capability::FittedBoundaryGeometry);
    if (request.ibm && !request.ibmForcing) {
        geometry.requiresAll.add(
            plan.kind == Kind::EulerianEulerian
                ? Capability::PhaseGhostState
                : Capability::ConservativeGhostState);
    } else if (plan.kind == Kind::EulerianEulerian) {
        geometry.provides.add(Capability::PhaseGhostState);
    }
    graph.add(std::move(geometry));

    ModuleDescriptor reconstruction;
    if (request.ilw) {
        reconstruction.name = "ILWBoundaryReconstruction";
        reconstruction.requiresAll = {
            Capability::CharacteristicEigenSystem,
            Capability::FiveVariableEulerState};
        reconstruction.provides = {
            Capability::ILWBoundaryReconstruction,
            Capability::ConservativeGhostState};
    } else {
        reconstruction.name = "AlgebraicBoundaryReconstruction";
        reconstruction.provides.add(
            Capability::AlgebraicBoundaryReconstruction);
        if (plan.kind != Kind::Homogeneous) {
            reconstruction.provides.add(
                Capability::ConservativeGhostState);
        }
        // 双欧拉的拟合边界由 PhaseBoundaryApplicator 闭合；当前 IBM 尚未提供
        // 每相 canonical ghost state，不能在能力图中把普通代数边界误报为该能力。
    }
    graph.add(std::move(reconstruction));

    if (request.ibmForcing) {
        ImmersedAlgorithm::validateForAlgorithm(
            config.ibm.forcing, config.numerics.solver);
        ModuleDescriptor constraint;
        const auto enforcement = config.ibm.forcing.enforcement;
        const ImmersedAlgorithm::StrategyContract strategy =
            ImmersedAlgorithm::contract(enforcement);
        constraint.name = strategy.name;
        constraint.requiresAll = {
            Capability::CanonicalConservativeState,
            Capability::ImmersedBoundaryGeometry};
        constraint.requiresAll.add(strategy.requiresPressureCorrection
            ? Capability::PressureCorrection
            : Capability::PredictedConservativeState);
        if (strategy.producesMultiplier) {
            constraint.provides.add(
                Capability::LagrangeMultiplierField);
        }
        if (strategy.requiresPredictedState
            || strategy.requiresPressureCorrection) {
            constraint.provides.add(Capability::ImmersedConstraintProjection);
        }
        graph.add(std::move(constraint));
    }

    if (request.interfaceGhostFluid) {
        ModuleDescriptor interface;
        interface.name = "GhostFluidInterface";
        interface.requiresOneOf.push_back({
            Capability::PressureJumpConsumer,
            Capability::TwoMaterialRiemann});
        graph.add(std::move(interface));
    }
    graph.validate();
    return graph;
}

} // namespace SF::Workflow
