/// @file SF_linearAlgebra.cpp
/// @brief 稀疏系统校验、结构复用会话与后端无关求解入口。

#include "solver/linearAlgebra/SF_linearAlgebra.h"

#include "solver/linearAlgebra/hypre/SF_hypre.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace SF::LinearAlgebra {

void validate(const SparseSystem& system) {
    if (system.globalSize <= 0 || system.firstRow < 0
        || system.lastRow < system.firstRow
        || system.lastRow >= system.globalSize) {
        throw std::runtime_error(
            "Linear sparse system has invalid global/local row bounds.");
    }
    const auto localSize = system.lastRow - system.firstRow + 1;
    if (static_cast<std::int64_t>(system.rows.size()) != localSize
        || static_cast<std::int64_t>(system.rhs.size()) != localSize
        || (!system.initialGuess.empty()
            && static_cast<std::int64_t>(system.initialGuess.size())
                != localSize)) {
        throw std::runtime_error(
            "Linear sparse system local arrays do not match owned rows.");
    }
    for (std::int64_t local = 0; local < localSize; ++local) {
        const SparseRow& row = system.rows[static_cast<size_t>(local)];
        if (row.globalRow != system.firstRow + local
            || row.columns.empty()
            || row.columns.size() != row.values.size()
            || !std::isfinite(system.rhs[static_cast<size_t>(local)])) {
            std::ostringstream message;
            message << "Linear sparse system contains an invalid local row: "
                    << "local=" << local << ", globalRow=" << row.globalRow
                    << ", expected=" << system.firstRow + local
                    << ", columns=" << row.columns.size()
                    << ", values=" << row.values.size()
                    << ", rhs="
                    << system.rhs[static_cast<size_t>(local)] << ".";
            throw std::runtime_error(message.str());
        }
        for (size_t n = 0; n < row.columns.size(); ++n) {
            if (row.columns[n] < 0 || row.columns[n] >= system.globalSize
                || !std::isfinite(row.values[n])) {
                throw std::runtime_error(
                    "Linear sparse system contains an invalid column/value.");
            }
        }
    }
}

SolverSession::SolverSession(FDM::LinearSolverConfig config)
    : config_(std::move(config)) {
    FDM::validateLinearSolverConfig(config, "linear solver");
    if (config_.backend != FDM::LinearBackend::Hypre) {
        throw std::runtime_error(
            "Requested linear backend is not implemented; use hypre.");
    }
    backend_ = std::make_unique<HypreSolverSession>(config_);
}

SolverSession::~SolverSession() = default;
SolverSession::SolverSession(SolverSession&&) noexcept = default;
SolverSession& SolverSession::operator=(SolverSession&&) noexcept = default;

SolveResult SolverSession::solve(const SparseSystem& system) {
    validate(system);
    if (config_.equilibration == FDM::LinearEquilibration::None) {
        return backend_->solve(system);
    }
    SparseSystem scaled=system;
    // 左行缩放只影响预条件和线性停止范数；不缩放解变量。
    for(size_t row=0;row<scaled.rows.size();++row) {
        double maximum=0.0;
        for(double value:scaled.rows[row].values)
            maximum=std::max(maximum,std::abs(value));
        if(maximum==0.0) {
            throw std::runtime_error(
                "rowMax equilibration found an all-zero matrix row.");
        }
        const double inverse=1.0/maximum;
        for(double& value:scaled.rows[row].values)value*=inverse;
        scaled.rhs[row]*=inverse;
    }
    validate(scaled);
    return backend_->solve(scaled);
}

void SolverSession::invalidateStructure() {
    backend_->invalidateStructure();
}

const ReuseStatistics& SolverSession::statistics() const {
    return backend_->statistics();
}

SolveResult solve(const SparseSystem& system,
                  const FDM::LinearSolverConfig& config) {
    SolverSession session(config);
    return session.solve(system);
}

} // namespace SF::LinearAlgebra
