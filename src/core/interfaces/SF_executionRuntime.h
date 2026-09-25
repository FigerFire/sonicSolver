#pragma once

/// @file SF_executionRuntime.h
/// @brief 求解算法可见的确定性执行契约；不暴露 MPI、rank 或通信计划。

#include "SF_stateBundle.h"

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace SF::Execution {

/// @brief 一个算子对分布式数据的访问语义。
enum class AccessMode {
    ReadOwned,             ///< 只读本地 owner 数据，不要求 halo。
    ReadHalo,              ///< 读取显式深度的最新 halo。
    WriteOwned,            ///< 写 owner 数据；完成后使该字段 halo 失效。
    AccumulateGlobalDof,   ///< 多个贡献累加到唯一 GlobalDof owner。
    WriteCanonicalFace     ///< 唯一 face owner 写通量，完成后装配 GlobalFace。
};

/// @brief 单个字段访问要求。
struct FieldAccess {
    std::string field;
    AccessMode mode = AccessMode::ReadOwned;
    int haloDepth = 0;

    void validate(const std::string& operation) const {
        if (field.empty()) {
            throw std::runtime_error(
                "Execution contract '" + operation
                + "' contains an empty field name.");
        }
        if (haloDepth < 0) {
            throw std::runtime_error(
                "Execution contract '" + operation + "' field '" + field
                + "' has a negative halo depth.");
        }
        if (mode == AccessMode::ReadHalo && haloDepth <= 0) {
            throw std::runtime_error(
                "Execution contract '" + operation + "' field '" + field
                + "' requires ReadHalo with an explicit positive depth.");
        }
        if (mode != AccessMode::ReadHalo && haloDepth != 0) {
            throw std::runtime_error(
                "Execution contract '" + operation + "' field '" + field
                + "' assigns haloDepth to a non-ReadHalo access.");
        }
    }
};

/// @brief 一个同步执行算子的完整数据契约。
struct OperatorContract {
    std::string name;
    std::vector<FieldAccess> accesses;

    void validate() const {
        if (name.empty()) {
            throw std::runtime_error(
                "Execution operator contract requires an explicit name.");
        }
        for (const auto& access : accesses) access.validate(name);
    }
};

inline FieldAccess readOwned(std::string field) {
    return {std::move(field), AccessMode::ReadOwned, 0};
}

inline FieldAccess readHalo(std::string field, int depth) {
    return {std::move(field), AccessMode::ReadHalo, depth};
}

inline FieldAccess writeOwned(std::string field) {
    return {std::move(field), AccessMode::WriteOwned, 0};
}

inline FieldAccess accumulateGlobalDof(std::string field) {
    return {std::move(field), AccessMode::AccumulateGlobalDof, 0};
}

inline FieldAccess writeCanonicalFace(std::string field) {
    return {std::move(field), AccessMode::WriteCanonicalFace, 0};
}

} // namespace SF::Execution

namespace SF { class FluxField; class Residual; }

namespace SF::FDM {

/// @brief 连续全局编号中当前执行单元拥有的区间。
struct DistributedIndexRange {
    std::int64_t first = 0;
    std::int64_t last = -1;
    std::int64_t total = 0;
};

/// @brief Algorithm 与数据分布实现之间的同步执行边界。
///
/// Algorithm 仍按数学顺序调用 `prepare -> operator -> finalize`。实现只根据显式
/// contract 和字段 freshness 决定是否交换 halo 或装配全局贡献；任何字段访问器
/// 都不会隐式触发通信。
class IExecutionRuntime {
public:
    virtual ~IExecutionRuntime() = default;
    virtual void attachState(State::StateBundle& state) = 0;
    virtual void prepare(const Execution::OperatorContract& contract) = 0;
    virtual void finalize(const Execution::OperatorContract& contract) = 0;
    virtual void synchronizeCanonicalFaceFluxes(
        const std::vector<Field*>& fields,
        const std::vector<FluxField*>& fluxes,
        const std::vector<Residual*>& residuals) = 0;
    virtual void accumulateGlobalDofResiduals(
        const std::vector<Field*>& fields,
        const std::vector<Residual*>& residuals) = 0;
    virtual double globalMinimum(double localValue) = 0;
    virtual double globalMaximum(double localValue) = 0;
    virtual double globalSum(double localValue) = 0;
    /// @brief 显式完成短生命周期数值向量的逐项 SUM，不注册为 Field。
    virtual void globalSum(std::vector<double>& values) = 0;
    /// @brief 查询一个全局实体是否由当前执行单元负责其唯一标量输出。
    ///
    /// 用于 marker 诊断、刚体载荷和 owner-only matrix row；它不暴露 rank，
    /// 更不允许物理模型自行组织通信。实体编号必须稳定且非负。
    virtual bool ownsCanonicalEntity(std::int64_t entityId) const = 0;
    /// @brief 当前 Runtime 是否连接到活动的分布式执行后端。
    virtual bool distributed() const = 0;
    /// @brief 将 owner 写入的 canonical 向量显式稀疏复制到所有 replicas。
    ///
    /// 调用者必须只在 `ownsCanonicalEntity(id)` 为真时写入非零值；Runtime
    /// 会在通信边界校验该契约并完成 owner-to-copy。该操作不触发 Field 访问器。
    virtual void copyCanonicalEntities(
        const std::vector<std::int64_t>& entityIds,
        std::vector<double>& values) = 0;
    /// @brief 为 owner 自由度分配连续、无重叠的分布式编号区间。
    virtual DistributedIndexRange allocateDistributedIndices(
        std::int64_t localCount) = 0;
    /// @brief 同步短生命周期算法 workspace，不把它注册为长期 Field。
    virtual void synchronizeTransient(
        const std::vector<State::DistributedFieldView>& fields) = 0;
    /// @brief 从重复 patch 候选中选择 canonical Eulerian owner cell。
    virtual std::vector<int> canonicalOwnerCells(
        Field& geometry,
        const std::vector<int>& candidates) = 0;
};


} // namespace SF::FDM
