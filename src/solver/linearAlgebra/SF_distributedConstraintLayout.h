#pragma once

/// @file SF_distributedConstraintLayout.h
/// @brief 将 ConstraintGlobalDof 映射到 owner-only 的连续线性代数行。

#include "core/interfaces/SF_executionRuntime.h"
#include "core/model/SF_globalEntityId.h"

#include <cstdint>
#include <unordered_map>
#include <vector>

namespace SF::LinearAlgebra {

/// @brief MPI/HYPRE 无关的约束行布局。
///
/// 所有 rank 可以持有相同的 marker identity 列表；只有 canonical owner
/// 分配本地连续行，随后通过 Runtime 的显式 owner-to-copy 让所有 replica
/// 获得同一 row。该类不保存 communicator，也不把 GlobalDof 与 backend row
/// 混为一类。
class DistributedConstraintLayout final {
public:
    /// @brief 根据稳定 ConstraintGlobalDof 建立一次不可变 row layout。
    static DistributedConstraintLayout build(
        const std::vector<GlobalConstraintDofId>& entities,
        FDM::IExecutionRuntime& runtime);

    /// @brief 查询任意已登记约束实体对应的 backend row。
    std::int64_t row(GlobalConstraintDofId entity) const;
    std::int64_t firstRow() const { return firstRow_; }
    std::int64_t lastRow() const { return lastRow_; }
    std::int64_t globalSize() const { return globalSize_; }
    const std::vector<GlobalConstraintDofId>& ownedEntities() const {
        return ownedEntities_;
    }

private:
    std::int64_t firstRow_ = 0;
    std::int64_t lastRow_ = -1;
    std::int64_t globalSize_ = 0;
    std::unordered_map<std::int64_t,std::int64_t> rows_;
    std::vector<GlobalConstraintDofId> ownedEntities_;
};

} // namespace SF::LinearAlgebra
