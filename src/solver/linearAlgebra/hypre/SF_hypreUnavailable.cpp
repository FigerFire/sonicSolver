/// @file SF_hypreUnavailable.cpp
/// @brief 未链接 HYPRE 时提供明确失败的后端占位实现。

#include "solver/linearAlgebra/hypre/SF_hypre.h"

#include <stdexcept>
#include <utility>

namespace SF::LinearAlgebra {

bool hypreBackendLinked() { return false; }

struct HypreSolverSession::Impl {
    explicit Impl(FDM::LinearSolverConfig value)
        : config(std::move(value)) {}

    FDM::LinearSolverConfig config;
    ReuseStatistics statistics;
};

HypreSolverSession::HypreSolverSession(FDM::LinearSolverConfig config)
    : impl_(std::make_unique<Impl>(std::move(config))) {
    throw std::runtime_error(
        "This workflow requires HYPRE with MPI, but sonicSolver was built "
        "without that optional backend. Install HYPRE/MPI and rebuild, or "
        "select a density-based workflow that does not assemble a global "
        "linear system.");
}

HypreSolverSession::~HypreSolverSession() = default;
HypreSolverSession::HypreSolverSession(HypreSolverSession&&) noexcept =
    default;
HypreSolverSession& HypreSolverSession::operator=(
    HypreSolverSession&&) noexcept = default;

SolveResult HypreSolverSession::solve(const SparseSystem&) {
    throw std::runtime_error(
        "HYPRE solve requested from a build without HYPRE/MPI.");
}

void HypreSolverSession::invalidateStructure() {}

const ReuseStatistics& HypreSolverSession::statistics() const {
    return impl_->statistics;
}

} // namespace SF::LinearAlgebra
