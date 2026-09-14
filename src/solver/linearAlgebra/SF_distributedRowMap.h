#pragma once

/// @file SF_distributedRowMap.h
/// @brief 压力基迁移期间的分布式矩阵行号适配层。

#include "core/field/SF_field.h"
#include "core/interfaces/SF_executionRuntime.h"

#include <cstdint>
#include <functional>
#include <vector>

namespace SF::LinearAlgebra {

/// @brief 当前 HYPRE SparseSystem 所需的本地 cell 到分布式 row 映射。
///
/// 这是压力基迁移期间的 legacy 适配类型，不是 GlobalDofId。新的方程装配不得
/// 读取 first/last/rank；目标接口是先按 GlobalDofId 形成 stencil，再由
/// DistributedLinearSystem 在 infrastructure 后端编号。
struct DistributedRowMap {
    std::vector<std::int64_t> global;
    std::vector<int> localCells;
    std::int64_t first = 0;
    std::int64_t last = -1;
    std::int64_t total = 0;
};

using CellSelector = std::function<bool(int, int, int)>;

/// @brief 为旧压力矩阵构造分布式行映射；ownership/通信由 Runtime 完成。
DistributedRowMap makeDistributedRowMap(
    const Field& field,
    FDM::IExecutionRuntime* runtime,
    const CellSelector& selected);

} // namespace SF::LinearAlgebra
