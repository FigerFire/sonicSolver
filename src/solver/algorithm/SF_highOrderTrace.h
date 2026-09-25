#pragma once

// 临时高阶回归追踪：只读现有 state/workspace，默认关闭（SF_HIGH_ORDER_TRACE=1）。

#include "core/field/SF_field.h"
#include "core/flux/SF_flux.h"
#include "core/residual/SF_residual.h"
#include "solver/algorithm/SF_patchWorkspace.h"
#include "methods/numerics/structured/SF_structured.h"
#include "core/state/SF_fluidStateModel.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace SF::SolverAlgorithm::HighOrderTrace {

inline bool enabled() {
    static const bool value = [] {
        const char* env = std::getenv("SF_HIGH_ORDER_TRACE");
        return env && std::string(env) == "1";
    }();
    return value;
}

inline std::int64_t selectedStep() {
    static const std::int64_t value = [] {
        const char* env = std::getenv("SF_HIGH_ORDER_TRACE_STEP");
        if (!env || *env == '\0') return std::int64_t{0};
        char* end = nullptr;
        errno = 0;
        const long long parsed = std::strtoll(env, &end, 10);
        if (errno != 0 || end == env || *end != '\0' || parsed < 0) {
            throw std::runtime_error(
                "SF_HIGH_ORDER_TRACE_STEP must be a non-negative integer.");
        }
        return static_cast<std::int64_t>(parsed);
    }();
    return value;
}

struct Context {
    std::int64_t step = -1;
    int stage = -1;
    double stageTime = 0.0;
};

inline Context& context() {
    static Context value;
    return value;
}

inline void beginPhysicalStep(std::int64_t step) {
    if (!enabled()) return;
    context() = {step, -1, 0.0};
}

inline void beginStage(std::int64_t step, double stageTime) {
    if (!enabled()) return;
    Context& value = context();
    if (value.step != step) value = {step, -1, 0.0};
    ++value.stage;
    value.stageTime = stageTime;
}

inline bool activeFor(std::int64_t step) {
    return enabled() && step == selectedStep();
}

inline std::ostream& printPrefix(std::ostream& output) {
    output << "[SF TRACE]";
    const Context& value = context();
    if (value.step >= 0) output << " step=" << value.step;
    if (value.stage >= 0) {
        output << " stage=" << value.stage
               << " stage-time=" << std::setprecision(17) << value.stageTime;
    }
    return output << ' ';
}

struct Statistics {
    std::size_t count = 0;
    double sum = 0.0;
    double l1 = 0.0;
    double l2 = 0.0;
    double minimum = std::numeric_limits<double>::infinity();
    double maximum = -std::numeric_limits<double>::infinity();
    double checksum = 0.0;

    void add(double value, std::size_t index) {
        ++count;
        sum += value;
        l1 += std::abs(value);
        l2 += value * value;
        minimum = std::min(minimum, value);
        maximum = std::max(maximum, value);
        checksum += static_cast<double>(index + 1) * value;
    }
};

inline void printStatistics(const char* label, int component,
                            const Statistics& values) {
    printPrefix(std::cout << std::setprecision(17))
              << label << " component=" << component
              << " count=" << values.count
              << " sum=" << values.sum
              << " l1=" << values.l1
              << " l2=" << std::sqrt(values.l2)
              << " min=" << values.minimum
              << " max=" << values.maximum
              << " checksum=" << values.checksum << '\n';
}

inline void conservative(const char* checkpoint,
                         const std::vector<Field*>& fields,
                         bool includeGhost = false) {
    if (!enabled()) return;
    for (std::size_t patch = 0; patch < fields.size(); ++patch) {
        const Field* field = fields[patch];
        if (!field) continue;
        std::vector<Statistics> q(static_cast<std::size_t>(field->NVar()));
        Statistics pressure;
        std::vector<Statistics> ghostQ(static_cast<std::size_t>(field->NVar()));
        Statistics ghostPressure;
        const auto add = [&](int i, int j, int k) {
            const std::size_t point = static_cast<std::size_t>(field->getIdx(i, j, k));
            for (int v = 0; v < field->NVar(); ++v) {
                q[static_cast<std::size_t>(v)].add((*field)(i, j, k, v), point);
            }
            pressure.add(field->thermodynamicState(i, j, k).pressure, point);
        };
        if (includeGhost) {
            for (int k = 0; k < field->MZ(); ++k)
                for (int j = 0; j < field->MY(); ++j)
                    for (int i = 0; i < field->MX(); ++i) {
                        add(i, j, k);
                        if (i < field->NG() || i >= field->MX() - field->NG()
                            || j < field->NG() || j >= field->MY() - field->NG()
                            || k < field->NG() || k >= field->MZ() - field->NG()) {
                            const std::size_t point = static_cast<std::size_t>(
                                field->getIdx(i, j, k));
                            for (int v = 0; v < field->NVar(); ++v) {
                                ghostQ[static_cast<std::size_t>(v)].add(
                                    (*field)(i, j, k, v), point);
                            }
                            ghostPressure.add(
                                field->thermodynamicState(i, j, k).pressure, point);
                        }
                    }
        } else {
            Math::forFluidInterior(*field, add);
        }
        printPrefix(std::cout) << checkpoint << " patch=" << patch
                  << " field=" << static_cast<const void*>(field)
                  << " region=" << (includeGhost ? "storage" : "interior") << '\n';
        for (int v = 0; v < field->NVar(); ++v) {
            printStatistics("Q", v, q[static_cast<std::size_t>(v)]);
        }
        printStatistics("pressure", 0, pressure);
        if (includeGhost) {
            for (int v = 0; v < field->NVar(); ++v) {
                printStatistics("ghost-Q", v, ghostQ[static_cast<std::size_t>(v)]);
            }
            printStatistics("ghost-pressure", 0, ghostPressure);
        }
        const int i0 = field->NG();
        const int i1 = field->MX() - field->NG() - 1;
        const int j = field->NG();
        const int k = field->NG();
        printPrefix(std::cout) << checkpoint << " patch=" << patch
                  << " x-samples physical-left=";
        for (int v = 0; v < field->NVar(); ++v) std::cout << (*field)(i0, j, k, v) << ' ';
        std::cout << " physical-right=";
        for (int v = 0; v < field->NVar(); ++v) std::cout << (*field)(i1, j, k, v) << ' ';
        std::cout << " ghost-left=";
        for (int v = 0; v < field->NVar(); ++v) std::cout << (*field)(i0 - 1, j, k, v) << ' ';
        std::cout << " ghost-right=";
        for (int v = 0; v < field->NVar(); ++v) std::cout << (*field)(i1 + 1, j, k, v) << ' ';
        std::cout << '\n';
    }
}

inline void flux(const char* checkpoint, const std::vector<Field*>& fields,
                 const std::vector<PatchWorkspace*>& workspaces) {
    if (!enabled()) return;
    for (std::size_t patch = 0; patch < fields.size(); ++patch) {
        const Field* field = fields[patch];
        const PatchWorkspace* workspace = workspaces[patch];
        if (!field || !workspace) continue;
        const FluxField& flux = workspace->convectiveFlux;
        printPrefix(std::cout) << checkpoint << " patch=" << patch
                  << " faces=" << flux.faceCount() << " nVar=" << flux.variableCount() << '\n';
        for (int dir = 0; dir < 3; ++dir) {
            for (int v = 0; v < flux.variableCount(); ++v) {
                Statistics values;
                for (int cell = 0; cell < field->TotalSize(); ++cell) {
                    const std::size_t face = static_cast<std::size_t>(dir)
                        * static_cast<std::size_t>(field->TotalSize())
                        + static_cast<std::size_t>(cell);
                    values.add(flux(face, v), static_cast<std::size_t>(cell));
                }
                printStatistics(dir == 0 ? "flux-XI" : dir == 1 ? "flux-ETA" : "flux-ZETA",
                                v, values);
            }
        }
        const int j = field->NG();
        const int k = field->NG();
        const int faces[] = {field->NG(), field->NG() + field->NX() / 2 - 1,
                             field->NG() + field->NX() / 2,
                             field->NG() + field->NX() / 2 + 1};
        for (int i : faces) {
            const std::size_t face = static_cast<std::size_t>(field->getIdx(i, j, k));
            printPrefix(std::cout) << checkpoint << " patch=" << patch
                      << " fixed-face=XI(" << i << ',' << j << ',' << k << ") F=";
            for (int v = 0; v < flux.variableCount(); ++v) std::cout << flux(face, v) << ' ';
            std::cout << '\n';
        }
    }
}

inline void residual(const char* checkpoint, const std::vector<Field*>& fields,
                     const std::vector<PatchWorkspace*>& workspaces,
                     bool local, bool global) {
    if (!enabled()) return;
    for (std::size_t patch = 0; patch < fields.size(); ++patch) {
        const Field* field = fields[patch];
        const PatchWorkspace* workspace = workspaces[patch];
        if (!field || !workspace) continue;
        const Residual& residual = workspace->residual;
        std::vector<Statistics> x(static_cast<std::size_t>(field->NVar()));
        std::vector<Statistics> y(static_cast<std::size_t>(field->NVar()));
        std::vector<Statistics> z(static_cast<std::size_t>(field->NVar()));
        std::vector<Statistics> localValues(static_cast<std::size_t>(field->NVar()));
        std::vector<Statistics> globalValues(static_cast<std::size_t>(field->NVar()));
        std::size_t globalMask = 0;
        Math::forFluidInterior(*field, [&](int i, int j, int k) {
            const std::size_t point = static_cast<std::size_t>(field->getIdx(i, j, k));
            for (int v = 0; v < field->NVar(); ++v) {
                x[static_cast<std::size_t>(v)].add(residual.x(i, j, k, v), point);
                y[static_cast<std::size_t>(v)].add(residual.y(i, j, k, v), point);
                z[static_cast<std::size_t>(v)].add(residual.z(i, j, k, v), point);
                if (local) localValues[static_cast<std::size_t>(v)].add(
                    residual.local(i, j, k, v), point);
                if (global && residual.hasGlobal(i, j, k)) {
                    globalValues[static_cast<std::size_t>(v)].add(
                        residual.global(i, j, k, v), point);
                }
            }
            if (global && residual.hasGlobal(i, j, k)) ++globalMask;
        });
        std::cout << "[SF TRACE] " << checkpoint << " patch=" << patch << '\n';
        for (int v = 0; v < field->NVar(); ++v) {
            printStatistics("residual-XI", v, x[static_cast<std::size_t>(v)]);
            printStatistics("residual-ETA", v, y[static_cast<std::size_t>(v)]);
            printStatistics("residual-ZETA", v, z[static_cast<std::size_t>(v)]);
            if (local) printStatistics("residual-local", v, localValues[static_cast<std::size_t>(v)]);
            if (global) printStatistics("residual-global", v, globalValues[static_cast<std::size_t>(v)]);
        }
        if (global) std::cout << "[SF TRACE] " << checkpoint << " patch=" << patch
                              << " globalMask=" << globalMask << '\n';
        const int j = field->NG();
        const int k = field->NG();
        const int cells[] = {field->NG(), field->NG() + field->NX() / 2 - 1,
                             field->NG() + field->NX() / 2,
                             field->NG() + field->NX() / 2 + 1};
        for (int i : cells) {
            printPrefix(std::cout) << checkpoint << " patch=" << patch
                      << " fixed-cell=(" << i << ',' << j << ',' << k << ")"
                      << " residual=";
            for (int v = 0; v < field->NVar(); ++v) {
                std::cout << '[' << residual.x(i, j, k, v) << ','
                          << residual.y(i, j, k, v) << ','
                          << residual.z(i, j, k, v) << ']';
            }
            std::cout << '\n';
        }
        if (global) {
            std::size_t reported = 0;
            Math::forFluidInterior(*field, [&](int i, int j, int k) {
                if (reported == 4 || !residual.hasGlobal(i, j, k)
                    || field->globalDofId(i, j, k) < 0) return;
                printPrefix(std::cout) << checkpoint << " patch=" << patch
                          << " GlobalDofId=" << field->globalDofId(i, j, k)
                          << " owner-rank=" << field->globalDofOwnerRank(i, j, k)
                          << " local=";
                for (int v = 0; v < field->NVar(); ++v) {
                    std::cout << residual.local(i, j, k, v) << ' ';
                }
                std::cout << " global=";
                for (int v = 0; v < field->NVar(); ++v) {
                    std::cout << residual.global(i, j, k, v) << ' ';
                }
                std::cout << '\n';
                ++reported;
            });
        }
    }
}

struct ZeroStatistics {
    std::size_t count = 0;
    std::size_t nonZeroCount = 0;
    double l1 = 0.0;
    double linf = 0.0;
    bool hasFirstNonZero = false;
    std::string firstLocation;
    double firstValue = 0.0;

    void add(double value, std::string location) {
        ++count;
        const double magnitude = std::abs(value);
        l1 += magnitude;
        linf = std::max(linf, magnitude);
        if (value == 0.0) return;
        ++nonZeroCount;
        if (!hasFirstNonZero) {
            hasFirstNonZero = true;
            firstLocation = std::move(location);
            firstValue = value;
        }
    }
};

inline void printZeroStatistics(const char* storage, const ZeroStatistics& values) {
    printPrefix(std::cout << std::setprecision(17))
        << "workspace-after-clear " << storage
        << " count=" << values.count
        << " l1=" << values.l1
        << " linf=" << values.linf
        << " non-zero-count=" << values.nonZeroCount;
    if (values.hasFirstNonZero) {
        std::cout << " first-non-zero=" << values.firstLocation
                  << " value=" << values.firstValue;
    }
    std::cout << '\n';
}

/// Flux/Residual are scratch storage.  A non-zero value immediately after
/// System::begin means a stale stage can enter the next numerical assembly.
inline void workspaceAfterClear(const Field& field,
                                const PatchWorkspace& workspace,
                                std::size_t patch) {
    if (!enabled()) return;

    const FluxField& flux = workspace.convectiveFlux;
    const Residual& residual = workspace.residual;
    ZeroStatistics fluxValues;
    ZeroStatistics residualValues;
    for (std::size_t face = 0; face < flux.faceCount(); ++face) {
        for (int variable = 0; variable < flux.variableCount(); ++variable) {
            fluxValues.add(flux(face, variable),
                "face=" + std::to_string(face)
                + ",component=" + std::to_string(variable));
        }
    }
    for (int k = 0; k < field.MZ(); ++k) {
        for (int j = 0; j < field.MY(); ++j) {
            for (int i = 0; i < field.MX(); ++i) {
                const std::string point = "cell=(" + std::to_string(i) + ','
                    + std::to_string(j) + ',' + std::to_string(k) + ')';
                if (residual.hasGlobal(i, j, k)) {
                    residualValues.add(1.0, point + ",global-mask");
                }
                for (int variable = 0; variable < field.NVar(); ++variable) {
                    const std::string component = ",component="
                        + std::to_string(variable);
                    residualValues.add(residual.x(i, j, k, variable),
                        point + ",XI" + component);
                    residualValues.add(residual.y(i, j, k, variable),
                        point + ",ETA" + component);
                    residualValues.add(residual.z(i, j, k, variable),
                        point + ",ZETA" + component);
                    residualValues.add(residual.source(i, j, k, variable),
                        point + ",source" + component);
                    residualValues.add(residual.local(i, j, k, variable),
                        point + ",local" + component);
                    if (residual.hasGlobal(i, j, k)) {
                        residualValues.add(residual.global(i, j, k, variable),
                            point + ",global" + component);
                    }
                }
            }
        }
    }
    printPrefix(std::cout) << "workspace-after-clear patch=" << patch
        << " field=" << static_cast<const void*>(&field)
        << " workspace=" << static_cast<const void*>(&workspace)
        << " flux=" << static_cast<const void*>(&flux)
        << " residual=" << static_cast<const void*>(&residual) << '\n';
    printZeroStatistics("FluxField", fluxValues);
    printZeroStatistics("Residual", residualValues);
    if (fluxValues.nonZeroCount != 0 || residualValues.nonZeroCount != 0) {
        std::ostringstream message;
        message << "Workspace clear invariant failed at step=" << context().step
                << ", stage=" << context().stage
                << ", patch=" << patch
                << ", field=" << static_cast<const void*>(&field)
                << ", workspace=" << static_cast<const void*>(&workspace)
                << ", flux=" << static_cast<const void*>(&flux)
                << ", residual=" << static_cast<const void*>(&residual)
                << "; FluxField nonZeroCount=" << fluxValues.nonZeroCount
                << ", Residual nonZeroCount=" << residualValues.nonZeroCount;
        if (fluxValues.hasFirstNonZero) {
            message << ", first FluxField entry=" << fluxValues.firstLocation
                    << " value=" << fluxValues.firstValue;
        }
        if (residualValues.hasFirstNonZero) {
            message << ", first Residual entry=" << residualValues.firstLocation
                    << " value=" << residualValues.firstValue;
        }
        throw std::runtime_error(message.str());
    }
}

inline void gamma(const FDM::SolverConfig& config,
                  const Physics::FluidStateModel::Model& equations) {
    if (!enabled()) return;
    const auto equationGamma = equations.perfectGasGamma();
    printPrefix(std::cout << std::setprecision(17))
              << "thermodynamics config-gamma=" << config.numerics.idealGasGamma
              << " equation-gamma=" << (equationGamma ? *equationGamma : -1.0)
              << " selected-divDispatch-gamma="
              << (equationGamma ? *equationGamma : config.numerics.idealGasGamma) << '\n';
}

} // namespace SF::SolverAlgorithm::HighOrderTrace
