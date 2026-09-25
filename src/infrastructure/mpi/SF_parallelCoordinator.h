#pragma once

/// @file SF_parallelCoordinator.h
/// @brief MPI/Halo 实现对 Solver 暴露的统一注册状态协调器。

#include "core/interfaces/SF_parallelCoordinatorContract.h"
#include "SF_domainDecomposition.h"
#include "SF_haloExchange.h"
#include "SF_processorFlux.h"
#include "SF_reduction.h"

namespace SF {
struct MeshBlockField;
namespace Parallel {

/// @brief 把 MPI、halo plan 和全局归约封装在 infrastructure 内。
class ParallelCoordinator final : public FDM::IParallelCoordinator {
public:
    ParallelCoordinator(
        HaloExchange* halo,
        std::vector<MeshBlockField>* blocks,
        const Backend::CommunicationBackend* backend,
        const ReductionService* reduction,
        const DomainDecomposition* decomposition,
        const ProcessorFlux* processorFlux,
        bool canonicalInterfaceFlux);

    void attachState(State::StateBundle& state) override;
    void synchronizeRegistered(State::HaloSyncStage stage) override;
    void ensureHalo(const std::string& name, int requiredDepth) override;
    void markModified(const std::string& name) override;
    void synchronizeTransient(
        const std::vector<State::DistributedFieldView>& fields) override;
    void assembleCanonicalInterfaceFluxes(const std::vector<Field*>& fields,
                                          const std::vector<FluxField*>& fluxes,
                                          const std::vector<Residual*>& residuals) override;
    void accumulateGlobalDof(const std::string& field,
        const std::vector<Field*>& fields,
        const std::vector<Residual*>& residuals) override;
    double reduce(double localValue, Reduction operation) override;
    void reduceSum(std::vector<double>& values) override;
    void copyCanonicalEntities(
        const std::vector<std::int64_t>& entityIds,
        std::vector<double>& values,
        const std::vector<unsigned char>& ownerMask) override;
    FDM::DistributedIndexRange distributeIndices(
        std::int64_t localCount) override;
    void broadcast(std::vector<double>& values, int root) override;
    bool allRanksAgree(bool localValue) override;
    bool active() const override;
    int rank() const override;
    int size() const override;
    void barrier() override;

    void configureBlocks(
        std::vector<MeshBlockField>* blocks,
        bool canonicalInterfaceFlux);

private:
    HaloExchange* halo_ = nullptr;
    std::vector<MeshBlockField>* blocks_ = nullptr;
    const Backend::CommunicationBackend* backend_ = nullptr;
    const ReductionService* reduction_ = nullptr;
    const DomainDecomposition* decomposition_ = nullptr;
    const ProcessorFlux* processorFlux_ = nullptr;
    State::StateBundle* state_ = nullptr;
    bool canonicalInterfaceFlux_ = false;

    void exchangeViews(
        const std::vector<State::DistributedFieldView*>& fields);
};

} // namespace Parallel
} // namespace SF
