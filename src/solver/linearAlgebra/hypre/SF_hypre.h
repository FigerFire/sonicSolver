#pragma once

/// @file SF_hypre.h
/// @brief HYPRE IJ/ParCSR 后端入口。

#include "solver/linearAlgebra/SF_linearAlgebra.h"

#include <memory>

namespace SF::LinearAlgebra {

/// @brief HYPRE IJ/ParCSR 会话；实现重建、更新、求解三阶段生命周期。
class HypreSolverSession {
public:
    explicit HypreSolverSession(FDM::LinearSolverConfig config);
    ~HypreSolverSession();
    HypreSolverSession(HypreSolverSession&&) noexcept;
    HypreSolverSession& operator=(HypreSolverSession&&) noexcept;
    HypreSolverSession(const HypreSolverSession&) = delete;
    HypreSolverSession& operator=(const HypreSolverSession&) = delete;

    SolveResult solve(const SparseSystem& system);
    void invalidateStructure();
    const ReuseStatistics& statistics() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace SF::LinearAlgebra
