/// @file SF_corrector.cpp
/// @brief 单相 SIMPLE/PISO/PIMPLE 压力—速度校正循环。

#include "solver/algorithm/pressureBased/SF_corrector.h"

#include "core/mesh/SF_dimension.h"
#include "methods/numerics/structured/SF_structured.h"
#include "solver/algorithm/pressureBased/SF_pressureJump.h"
#include "SF_thermodynamicClosure.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <sstream>
#include <stdexcept>

namespace SF::PressureBased {
namespace {

struct RowMap {
    std::vector<std::int64_t> global;
    std::vector<int> localCells;
    std::int64_t first = 0, last = -1, total = 0;
};

struct FieldRowMaps {
    std::vector<RowMap> fields;
    std::int64_t first = 0, last = -1, total = 0;
};

size_t index(const Field& field, int i, int j, int k) {
    return static_cast<size_t>(field.getIdx(i, j, k));
}

void offset(int axis, int sign, int& di, int& dj, int& dk) {
    di = dj = dk = 0;
    if (axis == 0) di = sign;
    else if (axis == 1) dj = sign;
    else dk = sign;
}

std::array<double, 3> metric(
        const Field& field, int axis, int i, int j, int k) {
    if (axis == 0) return {
        field.XiX(i,j,k), field.XiY(i,j,k), field.XiZ(i,j,k)};
    if (axis == 1) return {
        field.EtX(i,j,k), field.EtY(i,j,k), field.EtZ(i,j,k)};
    return {field.ZeX(i,j,k), field.ZeY(i,j,k), field.ZeZ(i,j,k)};
}

double spacing(const Field& field, int i, int j, int k,
               int ni, int nj, int nk) {
    const double dx = field.X(ni,nj,nk) - field.X(i,j,k);
    const double dy = field.Y(ni,nj,nk) - field.Y(i,j,k);
    const double dz = field.Z(ni,nj,nk) - field.Z(i,j,k);
    const double result = std::sqrt(dx*dx + dy*dy + dz*dz);
    if (!std::isfinite(result) || result <= 0.0) {
        throw std::runtime_error(
            "pressureBase found degenerate neighbour spacing.");
    }
    return result;
}

bool unknown(const Field& field, int i, int j, int k) {
    const int ng = field.NG();
    return i >= ng && i < ng + field.NX()
        && j >= ng && j < ng + field.NY()
        && k >= ng && k < ng + field.NZ()
        && field.CellFlag(i,j,k) == FLUID_CELL
        && !field.isSolverBoundaryPoint(i,j,k);
}

double density(const Field& field, int i, int j, int k) {
    const int count = field.hasEquationSet()
        ? field.equationSet()->densityVariableCount() : 1;
    double result = 0.0;
    for (int variable = 0; variable < count; ++variable) {
        result += field(i,j,k,variable);
    }
    if (!std::isfinite(result) || result <= 0.0) {
        throw std::runtime_error(
            "pressureBase found invalid density.");
    }
    return result;
}

std::array<int, 3> momentum(const Field& field) {
    if (field.hasEquationSet()) return {
        field.equationSet()->momentumIndex(0),
        field.equationSet()->momentumIndex(1),
        field.equationSet()->momentumIndex(2)};
    return {RU, RV, RW};
}

std::array<double, 3> velocity(
        const Field& field, const std::array<int,3>& momentumIndex,
        int i, int j, int k) {
    const double rho = density(field, i, j, k);
    return {field(i,j,k,momentumIndex[0]) / rho,
            field(i,j,k,momentumIndex[1]) / rho,
            field(i,j,k,momentumIndex[2]) / rho};
}

double divergence(const Field& field,
                  const std::array<int,3>& momentumIndex,
                  int i, int j, int k) {
    double result = 0.0;
    const auto center = velocity(field, momentumIndex, i, j, k);
    for (int axis = 0; axis < 3; ++axis) {
        if (!Math::isDirectionActiveIndex(axis)) continue;
        int ip = 0, jp = 0, kp = 0, im = 0, jm = 0, km = 0;
        offset(axis, 1, ip, jp, kp);
        offset(axis, -1, im, jm, km);
        const auto upper = field.CellFlag(i+ip,j+jp,k+kp) == SOLID_CELL
            ? center : velocity(field, momentumIndex, i+ip,j+jp,k+kp);
        const auto lower = field.CellFlag(i+im,j+jm,k+km) == SOLID_CELL
            ? center : velocity(field, momentumIndex, i+im,j+jm,k+km);
        const auto coordinateMetric = metric(field, axis, i, j, k);
        for (int component = 0; component < 3; ++component) {
            result += 0.5 * coordinateMetric[(size_t)component]
                    * (upper[(size_t)component]
                       - lower[(size_t)component]);
        }
    }
    return result;
}

FieldRowMaps makeRowMaps(
        const std::vector<Field*>& fields,
        FDM::IExecutionRuntime* runtime) {
    FieldRowMaps result;
    result.fields.resize(fields.size());
    std::int64_t local = 0;
    for (size_t block = 0; block < fields.size(); ++block) {
        if (!fields[block]) {
            throw std::runtime_error(
                "pressureBase received a null pressure field.");
        }
        const Field& field = *fields[block];
        auto& map = result.fields[block];
        map.global.assign((size_t)field.TotalSize(), -1);
        const int ng = field.NG();
        for (int k = ng; k < ng + field.NZ(); ++k) {
            for (int j = ng; j < ng + field.NY(); ++j) {
                for (int i = ng; i < ng + field.NX(); ++i) {
                    if (unknown(field, i, j, k)) {
                        map.localCells.push_back(field.getIdx(i,j,k));
                    }
                }
            }
        }
        local += (std::int64_t)map.localCells.size();
    }
    const FDM::DistributedIndexRange ownership = runtime
        ? runtime->allocateDistributedIndices(local)
        : FDM::DistributedIndexRange{0, local - 1, local};
    const std::int64_t offsetValue = ownership.first;
    const std::int64_t total = ownership.total;
    result.first = offsetValue;
    result.last = offsetValue + local - 1;
    result.total = total;
    std::vector<std::vector<double>> encoded(fields.size());
    std::int64_t next = offsetValue;
    for (size_t block = 0; block < fields.size(); ++block) {
        auto& map = result.fields[block];
        map.first = next;
        map.total = total;
        encoded[block].assign((size_t)fields[block]->TotalSize(), -1.0);
        for (int cell : map.localCells) {
            map.global[(size_t)cell] = next;
            encoded[block][(size_t)cell] = (double)next++;
        }
        map.last = next - 1;
    }
    if (runtime) {
        std::vector<State::DistributedFieldView> views;
        views.reserve(fields.size());
        for (size_t block = 0; block < fields.size(); ++block) {
            views.push_back(State::workspaceView(
                "pressureGlobalRow", (int)block, *fields[block],
                encoded[block], 1, fields[block]->NG(),
                State::HaloSyncStage::None,
                State::ExchangeKind::Identifier));
        }
        runtime->synchronizeTransient(views);
    }
    for (size_t block = 0; block < fields.size(); ++block) {
        for (size_t cell = 0; cell < encoded[block].size(); ++cell) {
            if (encoded[block][cell] >= 0.0) {
                const double rounded = std::round(encoded[block][cell]);
                if (std::abs(encoded[block][cell] - rounded) > 1.0e-9
                    || rounded > 9007199254740992.0) {
                    throw std::runtime_error(
                        "pressureBase received invalid remote row id.");
                }
                result.fields[block].global[cell] = (std::int64_t)rounded;
            }
        }
    }
    return result;
}

} // namespace

Corrector::Corrector(FDM::BoundaryConfig boundaries,
                     FDM::PressureCorrectionConfig config,
                     double idealGasGamma)
    : boundaries_(std::move(boundaries)), config_(std::move(config)),
      idealGasGamma_(idealGasGamma),
      linearSolver_(config_.workflow.pressure) {}

CorrectionSummary Corrector::correct(Field& field, double dt) {
    return correct(std::vector<Field*>{&field}, dt);
}

CorrectionSummary Corrector::correct(
        const std::vector<Field*>& fields, double dt) {
    if (!std::isfinite(dt) || dt <= 0.0) {
        throw std::runtime_error(
            "pressureBase requires finite positive dt.");
    }
    if (config_.workflow.type != FDM::SolverAlgorithm::PressureBased) {
        throw std::runtime_error(
            "pressureBase requires system/solverProperties.");
    }
    FDM::validateSolverProperties(config_.workflow);
    if (fields.size() > 1 && !runtime_) {
        throw std::runtime_error(
            "multi-field pressureBase requires Parallelism services.");
    }
    const FieldRowMaps maps = makeRowMaps(fields, runtime_);
    if (maps.last < maps.first || maps.total <= 0) {
        throw std::runtime_error(
            "pressureBase has no pressure unknowns.");
    }
    if (config_.workflow.referenceCell < 0
        || config_.workflow.referenceCell >= maps.total) {
        throw std::runtime_error(
            "pressureBase referenceCell is outside pressure matrix.");
    }

    std::vector<std::vector<unsigned char>> fixed(fields.size());
    for (size_t block = 0; block < fields.size(); ++block) {
        Field& field = *fields[block];
        fixed[block].assign((size_t)field.TotalSize(), 0);
        for (const auto& boundary : boundaries_.energyFromPressure) {
            if (boundary.type != FIXED_VALUE) continue;
            const auto set = field.getAllSets().find(boundary.name);
            if (set == field.getAllSets().end()) continue;
            for (int cell : set->second) {
                if (cell >= 0 && cell < field.TotalSize()) {
                    fixed[block][(size_t)cell] = 1;
                }
            }
        }
    }

    LinearAlgebra::SparseSystem matrix;
    matrix.firstRow = maps.first;
    matrix.lastRow = maps.last;
    matrix.globalSize = maps.total;
    const size_t localRows = (size_t)(maps.last - maps.first + 1);
    matrix.rows.reserve(localRows);
    matrix.rhs.reserve(localRows);
    matrix.initialGuess.assign(localRows, 0.0);
    CorrectionSummary summary;
    for (size_t block = 0; block < fields.size(); ++block) {
        Field& field = *fields[block];
        const RowMap& map = maps.fields[block];
        const auto momentumIndex = momentum(field);
        const FDM::IInterfaceJumpCondition* jumpCondition =
            interfaceJumpProvider_ ? interfaceJumpProvider_(field) : nullptr;
        for (int cell : map.localCells) {
            int i = 0, j = 0, k = 0;
            field.getIJK(cell, i, j, k);
            LinearAlgebra::SparseRow row;
            row.globalRow = map.global[(size_t)cell];
            const double div = divergence(field, momentumIndex, i, j, k);
            summary.maxDivergenceBefore = std::max(
                summary.maxDivergenceBefore, std::abs(div));
            ++summary.correctedCells;
            if (row.globalRow == config_.workflow.referenceCell) {
                row.columns = {row.globalRow};
                row.values = {1.0};
                matrix.rows.push_back(std::move(row));
                matrix.rhs.push_back(0.0);
                continue;
            }
            const double rho = density(field, i, j, k);
            const double soundSquared = field.hasEquationSet()
                ? std::pow(field.thermodynamicState(i,j,k).soundSpeed, 2)
                : idealGasGamma_ * Boundary::pressureAt(field,i,j,k) / rho;
            if (!std::isfinite(soundSquared) || soundSquared <= 0.0) {
                throw std::runtime_error(
                    "pressureBase found invalid compressibility.");
            }
            double diagonal = 1.0 / (rho * soundSquared * dt * dt);
            double rhs = -div / dt;
            for (int axis = 0; axis < 3; ++axis) {
                if (!Math::isDirectionActiveIndex(axis)) continue;
                for (int sign : {-1, 1}) {
                    int di = 0, dj = 0, dk = 0;
                    offset(axis, sign, di, dj, dk);
                    const int ni = i + di, nj = j + dj, nk = k + dk;
                    const int neighbour = field.getIdx(ni,nj,nk);
                    const auto global = map.global[(size_t)neighbour];
                    if (field.CellFlag(ni,nj,nk) != FLUID_CELL
                        || (global < 0
                            && !fixed[block][(size_t)neighbour])) continue;
                    const double h = spacing(field,i,j,k,ni,nj,nk);
                    const double coefficient = 0.5
                        * (1.0/rho + 1.0/density(field,ni,nj,nk))/(h*h);
                    const int leftI = sign > 0 ? i : ni;
                    const int leftJ = sign > 0 ? j : nj;
                    const int leftK = sign > 0 ? k : nk;
                    const FaceCorrectionJump jump = pressureCorrectionJump(
                        jumpCondition, field, leftI, leftJ, leftK, axis,
                        config_.workflow.pressureRelaxation);
                    if (jump.active) {
                        const double oriented = sign > 0
                            ? jump.correctionRightMinusLeft
                            : -jump.correctionRightMinusLeft;
                        // A p' 使用去跳跃差分：
                        // -c[(p'_n-p'_c)-J_cn]，故已知项移至 RHS 为 -c*J_cn。
                        rhs -= coefficient * oriented;
                        if (sign > 0) {
                            ++summary.interfaceFaces;
                            summary.maxTargetPressureJump = std::max(
                                summary.maxTargetPressureJump,
                                std::abs(jump.targetRightMinusLeft));
                            summary.maxCorrectionPressureJump = std::max(
                                summary.maxCorrectionPressureJump,
                                std::abs(jump.correctionRightMinusLeft));
                        }
                    }
                    diagonal += coefficient;
                    if (global >= 0) {
                        row.columns.push_back(global);
                        row.values.push_back(-coefficient);
                    }
                }
            }
            if (!std::isfinite(diagonal) || diagonal <= 0.0) {
                throw std::runtime_error(
                    "pressureBase produced invalid matrix diagonal.");
            }
            row.columns.push_back(row.globalRow);
            row.values.push_back(diagonal);
            matrix.rows.push_back(std::move(row));
            matrix.rhs.push_back(rhs);
            summary.initialResidual = std::max(
                summary.initialResidual, std::abs(rhs));
        }
    }

    const auto result = linearSolver_.solve(matrix);
    summary.iterations = result.iterations;
    summary.finalResidual = result.relativeResidual;
    summary.structureRebuilds =
        linearSolver_.statistics().structureRebuilds;
    summary.linearSolves = linearSolver_.statistics().solves;
    std::vector<std::vector<double>> correction(fields.size());
    size_t solutionRow = 0;
    for (size_t block = 0; block < fields.size(); ++block) {
        const RowMap& map = maps.fields[block];
        correction[block].assign((size_t)fields[block]->TotalSize(), 0.0);
        for (int cell : map.localCells) {
            correction[block][(size_t)cell] = result.solution[solutionRow++];
        }
    }
    if (runtime_) {
        std::vector<State::DistributedFieldView> views;
        views.reserve(fields.size());
        for (size_t block = 0; block < fields.size(); ++block) {
            views.push_back(State::workspaceView(
                "pressureCorrection", (int)block, *fields[block],
                correction[block], 1, fields[block]->NG(),
                State::HaloSyncStage::None));
        }
        runtime_->synchronizeTransient(views);
    }

    for (size_t block = 0; block < fields.size(); ++block) {
        Field& field = *fields[block];
        const RowMap& map = maps.fields[block];
        const auto momentumIndex = momentum(field);
        const FDM::IInterfaceJumpCondition* jumpCondition =
            interfaceJumpProvider_ ? interfaceJumpProvider_(field) : nullptr;
        for (int cell : map.localCells) {
            int i = 0, j = 0, k = 0;
            field.getIJK(cell, i, j, k);
            std::array<double,3> gradient{0.0,0.0,0.0};
            for (int axis = 0; axis < 3; ++axis) {
                if (!Math::isDirectionActiveIndex(axis)) continue;
                int ip=0,jp=0,kp=0,im=0,jm=0,km=0;
                offset(axis,1,ip,jp,kp); offset(axis,-1,im,jm,km);
                const auto coordinateMetric = metric(field,axis,i,j,k);
                const double upper = correction[block]
                    [index(field,i+ip,j+jp,k+kp)];
                const double lower = correction[block]
                    [index(field,i+im,j+jm,k+km)];
                double knownJump = 0.0;
                if (field.CellFlag(i+ip,j+jp,k+kp) == FLUID_CELL) {
                    const FaceCorrectionJump upperJump =
                        pressureCorrectionJump(
                            jumpCondition, field, i,j,k,axis,
                            config_.workflow.pressureRelaxation);
                    if (upperJump.active) {
                        knownJump += upperJump.correctionRightMinusLeft;
                    }
                }
                if (field.CellFlag(i+im,j+jm,k+km) == FLUID_CELL) {
                    const FaceCorrectionJump lowerJump =
                        pressureCorrectionJump(
                            jumpCondition, field,
                            i+im,j+jm,k+km,axis,
                            config_.workflow.pressureRelaxation);
                    if (lowerJump.active) {
                        knownJump += lowerJump.correctionRightMinusLeft;
                    }
                }
                for (int component = 0; component < 3; ++component) {
                    gradient[(size_t)component] += 0.5
                        * coordinateMetric[(size_t)component]
                        * (upper-lower-knownJump);
                }
            }
            for (int component = 0; component < 3; ++component) {
                field(i,j,k,momentumIndex[(size_t)component]) -=
                    config_.workflow.momentumRelaxation*dt
                    * gradient[(size_t)component];
            }
            const double pressure = Boundary::pressureAt(field,i,j,k)
                + config_.workflow.pressureRelaxation
                    * correction[block][(size_t)cell];
            if (!std::isfinite(pressure) || pressure <= 0.0) {
                throw std::runtime_error(
                    "pressureBase produced non-positive pressure.");
            }
            if (field.hasEquationSet()) {
                std::vector<double> state((size_t)field.NVar());
                for (int variable=0; variable<field.NVar(); ++variable)
                    state[(size_t)variable]=field(i,j,k,variable);
                field(i,j,k,field.equationSet()->energyIndex()) =
                    field.equationSet()->totalEnergyFromPressure(
                        state.data(),field.NVar(),pressure);
            } else {
                const double rho = density(field,i,j,k);
                double kinetic = 0.0;
                for (int component=0; component<3; ++component) {
                    const double value =
                        field(i,j,k,momentumIndex[(size_t)component]);
                    kinetic += 0.5*value*value/rho;
                }
                field(i,j,k,E)=pressure/(idealGasGamma_-1.0)+kinetic;
            }
        }
        for (int cell : map.localCells) {
            int i=0,j=0,k=0; field.getIJK(cell,i,j,k);
            summary.maxDivergenceAfter=std::max(
                summary.maxDivergenceAfter,
                std::abs(divergence(field,momentumIndex,i,j,k)));
        }
    }
    return summary;
}

} // namespace SF::PressureBased
