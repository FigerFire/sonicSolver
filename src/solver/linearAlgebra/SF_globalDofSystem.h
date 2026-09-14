#pragma once

/// @file SF_globalDofSystem.h
/// @brief MPI/HYPRE 无关的 GlobalDof 方程行与分布式线性系统边界。

#include "SF_linearAlgebra.h"

#include <cstdint>
#include <unordered_map>
#include <utility>
#include <vector>

namespace SF::LinearAlgebra {

/// @brief 数学自由度所属空间；空间标签防止不同物理未知量发生编号碰撞。
enum class GlobalDofSpace : std::uint8_t {
    Pressure = 1,
    Velocity = 2,
    Constraint = 3,
    Solid = 4,
    ScalarTransport = 5
};

/// @brief 网格/方程层的物理自由度身份；不是 HYPRE row id。
class GlobalDofId {
public:
    GlobalDofId() = default;
    explicit GlobalDofId(std::int64_t value) : value_(value) {}
    /// @brief 由未知量空间、canonical entity 和分量构造稳定数学身份。
    static GlobalDofId make(
        GlobalDofSpace space, std::int64_t entity, int component = 0);
    std::int64_t value() const { return value_; }
    bool valid() const { return value_ >= 0; }
    GlobalDofSpace space() const {
        return static_cast<GlobalDofSpace>((value_ >> 59) & 15);
    }
    friend bool operator==(GlobalDofId a, GlobalDofId b) {
        return a.value_ == b.value_;
    }
    friend bool operator!=(GlobalDofId a, GlobalDofId b) {
        return !(a == b);
    }
private:
    std::int64_t value_ = -1;
};

/// @brief 方程层的一项 `a(P,N)`，列身份仍是 GlobalDof。
struct GlobalDofCoefficient {
    GlobalDofId column;
    double value = 0.0;
};

/// @brief 方程层提交的一行，不包含 rank、ilower/iupper 或后端编号。
class GlobalDofRow {
public:
    explicit GlobalDofRow(GlobalDofId row = {}) : row_(row) {}
    GlobalDofId row() const { return row_; }
    const std::vector<GlobalDofCoefficient>& coefficients() const {
        return coefficients_;
    }
    void add(GlobalDofId column, double value);
    void setRightHandSide(double value) { rhs_ = value; }
    void setInitialGuess(double value) { initialGuess_ = value; }
    double rightHandSide() const { return rhs_; }
    double initialGuess() const { return initialGuess_; }
private:
    GlobalDofId row_;
    std::vector<GlobalDofCoefficient> coefficients_;
    double rhs_ = 0.0;
    double initialGuess_ = 0.0;
};

/// @brief 当前 rank 对 owner GlobalDof 装配的方程集合。
struct GlobalDofSystem {
    std::vector<GlobalDofRow> rows;
};

/// @brief Infrastructure 提供的 GlobalDof -> distributed row 映射。
class IDistributedNumbering {
public:
    virtual ~IDistributedNumbering() = default;
    virtual std::int64_t row(GlobalDofId dof) const = 0;
    virtual std::int64_t firstRow() const = 0;
    virtual std::int64_t lastRow() const = 0;
    virtual std::int64_t globalSize() const = 0;
};

/// @brief 已知 ownership 的不可变编号表；可由 MPI/HYPRE infrastructure 构造。
class StaticDistributedNumbering final : public IDistributedNumbering {
public:
    StaticDistributedNumbering(
        std::int64_t firstRow,
        std::int64_t lastRow,
        std::int64_t globalSize,
        const std::vector<std::pair<GlobalDofId, std::int64_t>>& entries);
    std::int64_t row(GlobalDofId dof) const override;
    std::int64_t firstRow() const override { return firstRow_; }
    std::int64_t lastRow() const override { return lastRow_; }
    std::int64_t globalSize() const override { return globalSize_; }
private:
    std::int64_t firstRow_ = 0;
    std::int64_t lastRow_ = -1;
    std::int64_t globalSize_ = 0;
    std::unordered_map<std::int64_t, std::int64_t> rows_;
};

/// @brief 唯一允许把 GlobalDof 行翻译成后端 SparseSystem 的层。
class DistributedLinearSystem {
public:
    DistributedLinearSystem(
        FDM::LinearSolverConfig config,
        const IDistributedNumbering& numbering);
    SolveResult solve(const GlobalDofSystem& equations);
    void invalidateStructure();
    const ReuseStatistics& statistics() const;
private:
    const IDistributedNumbering& numbering_;
    SolverSession session_;
};

} // namespace SF::LinearAlgebra
