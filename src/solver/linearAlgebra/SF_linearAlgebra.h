#pragma once

/// @file SF_linearAlgebra.h
/// @brief 分布式稀疏线性系统与求解服务接口。

#include "SF_config.h"

#include <cstdint>
#include <memory>
#include <vector>

namespace SF::LinearAlgebra {

/// @brief 后端已编号的一行稀疏矩阵数据。
/// @details 仅允许由 DistributedLinearSystem 或 legacy 适配器构造；这里的
/// globalRow/columns 是 distributed backend row，不是 GlobalDofId。
struct SparseRow {
    std::int64_t globalRow = -1;
    std::vector<std::int64_t> columns;
    std::vector<double> values;
};

/// @brief 一个 MPI rank 拥有的连续行分布线性系统。
struct SparseSystem {
    std::int64_t firstRow = 0;
    std::int64_t lastRow = -1;
    std::int64_t globalSize = 0;
    std::vector<SparseRow> rows;
    std::vector<double> rhs;
    std::vector<double> initialGuess;
    /// @brief 方程块语义：0=其他，1=压力，2=速度，3=约束，4=固体。
    std::vector<int> blockRoles;
};

/// @brief Krylov 求解聚合结果。
struct SolveResult {
    int iterations = 0;
    double relativeResidual = 0.0;
    bool converged = false;
    std::vector<double> solution;
};

/// @brief HYPRE 结构复用统计，用于确认重建/更新/求解的实际次数。
struct ReuseStatistics {
    std::int64_t structureRebuilds = 0;
    std::int64_t coefficientUpdates = 0;
    std::int64_t rhsUpdates = 0;
    std::int64_t preconditionerRefreshes = 0;
    std::int64_t solves = 0;
};

class HypreSolverSession;

/// @brief 可跨时间步和 corrector 复用矩阵拓扑及 HYPRE 求解器对象。
class SolverSession {
public:
    explicit SolverSession(FDM::LinearSolverConfig config);
    ~SolverSession();
    SolverSession(SolverSession&&) noexcept;
    SolverSession& operator=(SolverSession&&) noexcept;
    SolverSession(const SolverSession&) = delete;
    SolverSession& operator=(const SolverSession&) = delete;

    /// @brief 拓扑不变时仅更新 coefficients、rhs 和 initial guess。
    SolveResult solve(const SparseSystem& system);
    /// @brief 请求下一次 solve 重建矩阵结构和 AMG 层次。
    void invalidateStructure();
    const ReuseStatistics& statistics() const;

private:
    FDM::LinearSolverConfig config_;
    std::unique_ptr<HypreSolverSession> backend_;
};

/// @brief 校验本地稀疏系统的维数、所有权和有限性。
void validate(const SparseSystem& system);

/// @brief 使用配置声明的后端求解分布式线性系统。
SolveResult solve(const SparseSystem& system,
                  const FDM::LinearSolverConfig& config);

} // namespace SF::LinearAlgebra
