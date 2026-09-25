#pragma once

/// @file SF_communicationBackend.h
/// @brief infrastructure 并行模块使用的通信后端契约。

#include <cstdint>
#include <vector>

namespace SF::Parallel::Backend {

/// @brief 局部通信和全局通信的最小后端接口。
///
/// Field、Mesh 和 Solver 均不依赖该接口的具体实现。MPI communicator 以及
/// request/status 等对象只能出现在 backend 实现中。
class CommunicationBackend {
public:
    virtual ~CommunicationBackend() = default;

    virtual bool available() const = 0;
    virtual bool active() const = 0;
    virtual int rank() const = 0;
    virtual int size() const = 0;

    virtual bool allRanksAgree(bool localValue) const = 0;
    virtual double allReduceMin(double localValue) const = 0;
    virtual double allReduceMax(double localValue) const = 0;
    virtual double allReduceSum(double localValue) const = 0;
    /// @brief 对同一长度的数值向量执行逐项 SUM；只用于显式 Runtime 边界。
    virtual void allReduceSum(std::vector<double>& values) const = 0;
    /// @brief 按稳定 entity id 将 owner 值稀疏复制到所有 replicas。
    /// @param ownerMask 由上层 ownership registry 给出的唯一 owner 标记。
    virtual void copyCanonical(
        const std::vector<std::int64_t>& entityIds,
        std::vector<double>& values,
        const std::vector<unsigned char>& ownerMask) const = 0;
    virtual std::int64_t exclusiveScanSum(std::int64_t localValue) const = 0;
    virtual std::int64_t allReduceSum(std::int64_t localValue) const = 0;
    virtual void broadcast(std::vector<double>& values, int root) const = 0;
    virtual void barrier() const = 0;

    /// @brief 与给定邻居交换静态布局的双精度 payload。
    virtual void exchangeNeighbours(
        const std::vector<double>& send,
        const std::vector<int>& counts,
        const std::vector<int>& displacements,
        const std::vector<int>& neighbours,
        int tag,
        std::vector<double>& received) const = 0;
};

} // namespace SF::Parallel::Backend
