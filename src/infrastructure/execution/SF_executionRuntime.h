#pragma once

/// @file SF_executionRuntime.h
/// @brief 基于 freshness 和显式 contract 的同步执行 Runtime。

#include "core/interfaces/SF_executionRuntime.h"
#include "core/interfaces/SF_parallelCoordinatorContract.h"

namespace SF::Execution {

/// @brief 可用于单元测试和性能诊断的确定性 Runtime 计数器。
struct RuntimeStatistics {
    unsigned long long preparedOperators = 0;
    unsigned long long finalizedOperators = 0;
    unsigned long long haloRequests = 0;
    unsigned long long haloExchanges = 0;
    unsigned long long freshHaloReuses = 0;
    unsigned long long ownedWrites = 0;
    unsigned long long canonicalFaceAssemblies = 0;
    unsigned long long globalDofAssemblies = 0;
};

/// @brief 把抽象分布式服务隐藏在 Runtime 后方。
class Runtime final : public FDM::IExecutionRuntime {
public:
    explicit Runtime(FDM::IParallelCoordinator* parallel = nullptr)
        : parallel_(parallel) {}

    void attachState(State::StateBundle& state) override;
    void prepare(const OperatorContract& contract) override;
    void finalize(const OperatorContract& contract) override;
    void synchronizeCanonicalFaceFluxes(
        const std::vector<Field*>& fields,
        const std::vector<FluxField*>& fluxes,
        const std::vector<Residual*>& residuals) override;
    void accumulateGlobalDofResiduals(
        const std::vector<Field*>& fields,
        const std::vector<Residual*>& residuals) override;
    double globalMinimum(double localValue) override;
    double globalMaximum(double localValue) override;
    double globalSum(double localValue) override;
    void globalSum(std::vector<double>& values) override;
    bool ownsCanonicalEntity(std::int64_t entityId) const override;
    bool distributed() const override;
    void copyCanonicalEntities(
        const std::vector<std::int64_t>& entityIds,
        std::vector<double>& values) override;
    FDM::DistributedIndexRange allocateDistributedIndices(
        std::int64_t localCount) override;
    void synchronizeTransient(
        const std::vector<State::DistributedFieldView>& fields) override;
    std::vector<int> canonicalOwnerCells(
        Field& geometry,
        const std::vector<int>& candidates) override;

    const RuntimeStatistics& statistics() const { return statistics_; }

private:
    FDM::IParallelCoordinator* parallel_ = nullptr;
    State::StateBundle* state_ = nullptr;
    RuntimeStatistics statistics_;
};

} // namespace SF::Execution
