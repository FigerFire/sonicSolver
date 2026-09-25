#pragma once

/// @file SF_parallelCoordinatorContract.h
/// @brief Solver-visible ownership, synchronization, and reduction port.

#include "core/interfaces/SF_executionRuntime.h"

#include <cstdint>
#include <string>
#include <vector>

namespace SF::FDM {

class IParallelCoordinator {
public:
    virtual ~IParallelCoordinator() = default;
    virtual void attachState(State::StateBundle& state) = 0;
    virtual void synchronizeRegistered(State::HaloSyncStage stage) = 0;
    virtual void ensureHalo(const std::string& name, int requiredDepth) = 0;
    virtual void markModified(const std::string& name) = 0;
    virtual void synchronizeTransient(
        const std::vector<State::DistributedFieldView>& fields) = 0;
    virtual void assembleCanonicalInterfaceFluxes(
        const std::vector<Field*>& fields,
        const std::vector<FluxField*>& fluxes,
        const std::vector<Residual*>& residuals) = 0;
    virtual void accumulateGlobalDof(
        const std::string& field,
        const std::vector<Field*>& fields,
        const std::vector<Residual*>& residuals) = 0;

    enum class Reduction { Minimum, Maximum, Sum };
    virtual double reduce(double localValue, Reduction operation) = 0;
    virtual void reduceSum(std::vector<double>& values) = 0;
    virtual void copyCanonicalEntities(
        const std::vector<std::int64_t>& entityIds,
        std::vector<double>& values,
        const std::vector<unsigned char>& ownerMask) = 0;
    virtual DistributedIndexRange distributeIndices(
        std::int64_t localCount) = 0;
    virtual void broadcast(std::vector<double>& values, int root) = 0;
    virtual bool allRanksAgree(bool localValue) = 0;
    virtual bool active() const = 0;
    virtual int rank() const = 0;
    virtual int size() const = 0;
    virtual void barrier() = 0;
};

} // namespace SF::FDM
