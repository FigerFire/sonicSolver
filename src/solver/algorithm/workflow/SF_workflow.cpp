/// @file SF_workflow.cpp
/// @brief 从求解器、相系统和界面选择构造顶层工作流。

#include "SF_workflow.h"

#include <cmath>
#include <stdexcept>

namespace SF::Workflow {
namespace {

bool hasTransportedTurbulence(const FDM::TurbulenceConfig& turbulence) {
    return turbulence.enabled
        && turbulence.family != FDM::TurbulenceFamily::None
        && turbulence.family != FDM::TurbulenceFamily::DNS;
}

} // namespace

std::string Plan::name() const {
    switch (kind) {
        case Kind::SingleFluid:
            return "singleFluid";
        case Kind::OneFluidInterface:
            return "oneFluidInterface";
        case Kind::LegacyMultiphase:
            return "legacyMultiphase";
        case Kind::Homogeneous:
            return "homogeneous";
        case Kind::EulerianEulerian:
            return "eulerianEulerian";
    }
    throw std::runtime_error("Workflow plan contains an invalid kind.");
}

Plan makePlan(
        const Request& request,
        const FDM::SolverConfig& config) {
    if (request.homogeneous && request.eulerianEulerian) {
        throw std::runtime_error(
            "A case cannot select homogeneous and eulerianEulerian "
            "workflows simultaneously.");
    }

    Plan plan;
    plan.transportedLegacyAlpha=request.transportedLegacyAlpha;
    plan.turbulenceHasTransportState=
        hasTransportedTurbulence(config.turbulence);
    if (!request.multiphaseEnabled) {
        plan.kind=Kind::SingleFluid;
        return plan;
    }
    if (request.homogeneous) {
        plan.kind=Kind::Homogeneous;
        if (config.turbulence.enabled
            && config.turbulence.family != FDM::TurbulenceFamily::DNS) {
            throw std::runtime_error(
                "homogeneous EquationSet does not yet provide an "
                "EquationSet-aware turbulence system.");
        }
        if (config.numerics.solver
            != FDM::SolverAlgorithm::DensityBased) {
            throw std::runtime_error(
                "homogeneous EquationSet currently requires densityBased "
                "time integration.");
        }
        return plan;
    }
    if (request.eulerianEulerian) {
        plan.kind=Kind::EulerianEulerian;
        if (config.pressure.workflow.type
            != FDM::SolverAlgorithm::PressureBased) {
            throw std::runtime_error(
                "eulerianEulerian requires pressureBased workflow.");
        }
        if (!std::isfinite(config.numerics.maxDeltaT)
            || config.numerics.maxDeltaT<=0.0) {
            throw std::runtime_error(
                "eulerianEulerian requires finite positive maxDeltaT.");
        }
        return plan;
    }

    if (request.interfaceModel) {
        plan.kind = Kind::OneFluidInterface;
        if (request.phaseChange) {
            throw std::runtime_error(
                "OneFluid interface phase change is not implemented; "
                "select a supported homogeneous or Eulerian-Eulerian "
                "phase-change workflow.");
        }
        return plan;
    }

    plan.kind=Kind::LegacyMultiphase;
    if (request.phaseChange) {
        throw std::runtime_error(
            "Legacy multiphase phase change is not conservative. Select "
            "homogeneousMultiphase or eulerianEulerian explicitly.");
    }
    return plan;
}

void validateMeshExecution(
        const Plan& plan,
        bool multiField,
        bool hasCoupledInterfaces) {
    if (!multiField) return;
    if (plan.kind==Kind::Homogeneous) {
        throw std::runtime_error(
            "homogeneous EquationSet has not migrated to the multi-field "
            "MPI stepper.");
    }
    if (plan.kind==Kind::EulerianEulerian) {
        throw std::runtime_error(
            "eulerianEulerian requires one structured Field per MPI rank; "
            "multi-patch-per-rank execution is not implemented.");
    }
    if (plan.turbulenceHasTransportState) {
        throw std::runtime_error(
            "Multi-field MPI has not migrated transported turbulence scalar "
            "state and halo plans.");
    }
    if (plan.kind==Kind::LegacyMultiphase
        &&plan.transportedLegacyAlpha&&hasCoupledInterfaces) {
        throw std::runtime_error(
            "Legacy transported alpha on coupled patches has no canonical "
            "scalar interface flux.");
    }
}

} // namespace SF::Workflow
