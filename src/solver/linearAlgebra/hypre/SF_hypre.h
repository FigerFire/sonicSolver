#pragma once

/// @file SF_hypre.h
/// @brief HYPRE IJ/ParCSR 后端入口。

#include "solver/linearAlgebra/SF_linearAlgebra.h"

#include <memory>

namespace SF::LinearAlgebra {

/// @brief 本 binary 是否链接了真实 HYPRE ParCSR 后端。
///
/// 由两个互斥实现回答（真实后端 / 未链接占位实现），因此它反映的是
/// 二进制组成，而不是 case 或编译期假设。
bool hypreBackendLinked();

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
