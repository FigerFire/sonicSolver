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
        const std::string& flowLabel) {
    std::ostringstream os;
    os << config.type;
    if (!flowLabel.empty()) os << ", flow=" << flowLabel;
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
    broadcast("Flow formulation   : ",
              std::string("equations + coupling preset")
                  + (eulerianEulerian ? " (Eulerian shared pressure)"
                                      : ""));
    {
        const auto& pressure = config.pressure.linear.pressure;
        broadcast("Coupling algorithm : ",
                  FDM::toString(config.pressure.coupling.preset));
        broadcast("Pressure solve     : ",
                  "HYPRE maxIter=" + std::to_string(pressure.maxIterations)
                  + ", relTol=" + formatTimeValue(pressure.relativeTolerance)
                  + ", pRelax=" + formatTimeValue(
                      config.pressure.coupling.pressureRelaxation)
                  + ", uRelax=" + formatTimeValue(
                      config.pressure.coupling.momentumRelaxation));
    }
    broadcast("CFL number       : ", config.numerics.cfl);
    broadcast("Maximum deltaT   : ", config.numerics.maxDeltaT);
    broadcast("Active directions : ", Math::activeDirectionText());
    broadcast("Formulation      : ", FDM::toString(config.numerics.formulation));
    if (eulerianEulerian) {
        broadcast("Phase convection   : ",
                  FDM::toString(config.pressure.phaseTransport.convection));
        broadcast("Phase source CFL   : ",
                  config.pressure.phaseTransport.sourceCfl);
        broadcast("Time scheme       : ", "implicitEuler");
    } else {
        broadcast("Convection recipe : ",
                  config.numerics.recipes.convection
                    ? FDM::toString(config.numerics.recipes.convection->id())
                    : "unbound");
        broadcast("Interface flux    : ",
                  FDM::toString(config.numerics.interfaceFlux));
        broadcast("Time recipe       : ",
                  FDM::toString(config.numerics.timeRecipe.id()));
    }
    if (config.ibm.enabled
        && config.ibm.method == FDM::IBMMethod::Ghost) {
        broadcast("IBM WENO closure  : ",
                  FDM::toString(config.numerics.ibmBoundary));
    } else if (config.ibm.enabled) {
        broadcast("IBM WENO closure  : ", "not-applicable (forcing IBM)");
    }
    broadcast("ILW order         : ", config.numerics.ilwOrder);
    broadcast("Diffusion recipe  : ",
              config.numerics.recipes.diffusion
                ? FDM::toString(config.numerics.recipes.diffusion->id())
                : "not selected");
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
