/// @file SF_report.cpp
/// @brief 生成 solver configuration、run control 与 timestep 的日志文本。
///
/// Data flow:
///   resolved configuration + scalar runtime values
///       -> formatted strings / broadcast messages
///       -> CLI and solver log
///
/// 本文件不写 VTK、不改变求解状态，也不参与 execution routing。

#include "app/application/output/SF_report.h"

#include "core/interfaces/SF_log.h"
#include "core/mesh/SF_dimension.h"
#include "SF_multiphase.h"

#include <iomanip>
#include <iostream>
#include <sstream>

namespace SF::Application::Report {

std::string multiPhaseSummary(
        const Physics::Multiphase::MultiPhaseConfig& config,
        FDM::SolverAlgorithm solver) {
    std::ostringstream os;
    os << config.type << ", solver=" << FDM::toString(solver);
    if (!config.conservativeLayout.variables.empty()) {
        os << ", conserved=[";
        for (size_t i = 0; i < config.conservativeLayout.variables.size(); ++i) {
            if (i != 0) os << ",";
            os << config.conservativeLayout.variables[i];
        }
        os << "]";
    }
    return os.str();
}

std::string formatTimeValue(double value) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(6) << value;
    return oss.str();
}

std::string formatTimeStepStatus(double time, double dt) {
    return "time=" + formatTimeValue(time)
        + ", dt=" + formatTimeValue(dt);
}

void broadcastSolverConfig(
        const FDM::SolverConfig& config,
        bool eulerianEulerian) {
    if (logFromThisProcess()) std::cout << std::endl;
    broadcast("Solver             : ", FDM::toString(config.numerics.solver));
    if (config.numerics.solver == FDM::SolverAlgorithm::PressureBased) {
        const auto& pressure = config.pressure.workflow.pressure;
        broadcast("Coupling algorithm : ",
                  FDM::toString(config.pressure.workflow.algorithm));
        broadcast("Pressure solve     : ",
                  "HYPRE maxIter=" + std::to_string(pressure.maxIterations)
                  + ", relTol=" + formatTimeValue(pressure.relativeTolerance)
                  + ", pRelax=" + formatTimeValue(
                      config.pressure.workflow.pressureRelaxation)
                  + ", uRelax=" + formatTimeValue(
                      config.pressure.workflow.momentumRelaxation));
    }
    broadcast("CFL number       : ", config.numerics.cfl);
    broadcast("Maximum deltaT   : ", config.numerics.maxDeltaT);
    broadcast("Active directions : ", Math::activeDirectionText());
    broadcast("Formulation      : ", FDM::toString(config.numerics.formulation));
    if (eulerianEulerian) {
        broadcast("Phase convection   : ",
                  FDM::toString(config.pressure.workflow.phaseConvection));
        broadcast("Phase source CFL   : ",
                  config.pressure.workflow.phaseSourceCfl);
        broadcast("Time scheme       : ", "implicitEuler");
    } else {
        broadcast("Convection scheme : ",
                  FDM::toString(config.numerics.convection));
        broadcast("Reconstruction   : ",
                  FDM::toString(config.numerics.reconstruction));
        broadcast("Flux method       : ", FDM::toString(config.numerics.flux));
        broadcast("Interface flux    : ",
                  FDM::toString(config.numerics.interfaceFlux));
        broadcast("Time scheme       : ", FDM::toString(config.numerics.time));
    }
    if (config.ibm.enabled
        && config.ibm.method == FDM::IBMMethod::Ghost) {
        broadcast("IBM WENO closure  : ",
                  FDM::toString(config.numerics.ibmBoundary));
    } else if (config.ibm.enabled) {
        broadcast("IBM WENO closure  : ", "not-applicable (forcing IBM)");
    }
    broadcast("ILW order         : ", config.numerics.ilwOrder);
    broadcast("Viscous order     : ", FDM::toString(config.numerics.viscous));
    broadcast("Source terms      : ", FDM::toString(config.sources.enabled));
    if (config.turbulence.enabled) {
        broadcast("Turbulence family : ",
                  FDM::toString(config.turbulence.family));
        broadcast("Turbulence model  : ",
                  FDM::toString(config.turbulence.model));
    } else {
        broadcast("Turbulence model  : ", "Disabled");
    }
}

std::string formatRunControl(const CaseConfig& config) {
    const auto& time = config.time;
    std::ostringstream oss;
    oss << "startTime=" << formatTimeValue(time.startTime)
        << ", endTime=" << formatTimeValue(time.endTime)
        << ", endStep=" << time.endStep
        << ", output=" << (time.writeByStep ? "step" : "time")
        << ", saveBy=" << (time.writeByStep
                            ? std::to_string(time.writeIntervalSteps)
                            : formatTimeValue(time.writeIntervalTime))
        << ", writeInitial="
        << (config.writeInitial ? "true" : "false");
    return oss.str();
}

} // namespace SF::Application::Report
