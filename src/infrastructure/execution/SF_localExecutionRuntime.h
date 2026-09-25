#pragma once

/// @file SF_localExecutionRuntime.h
/// @brief Deterministic serial implementation of the execution contract.

#include "core/interfaces/SF_executionRuntime.h"

#include <cmath>
#include <stdexcept>

namespace SF::Execution {

class LocalRuntime final : public FDM::IExecutionRuntime {
public:
    void attachState(State::StateBundle& state) override { state_ = &state; }
    void prepare(const OperatorContract& contract) override;
    void finalize(const OperatorContract& contract) override;
    void synchronizeCanonicalFaceFluxes(
        const std::vector<Field*>&,
        const std::vector<FluxField*>&,
        const std::vector<Residual*>&) override {}
    void accumulateGlobalDofResiduals(
        const std::vector<Field*>&,
        const std::vector<Residual*>&) override {}
    double globalMinimum(double value) override { return value; }
    double globalMaximum(double value) override { return value; }
    double globalSum(double value) override { return value; }
    void globalSum(std::vector<double>&) override {}
    bool ownsCanonicalEntity(std::int64_t entityId) const override;
    bool distributed() const override { return false; }
    void copyCanonicalEntities(
        const std::vector<std::int64_t>& entityIds,
        std::vector<double>& values) override;
    FDM::DistributedIndexRange allocateDistributedIndices(
        std::int64_t localCount) override;
    void synchronizeTransient(
        const std::vector<State::DistributedFieldView>&) override {}
    std::vector<int> canonicalOwnerCells(
        Field&, const std::vector<int>& candidates) override {
        return candidates;
    }

private:
    State::StateBundle* state_ = nullptr;
};

} // namespace SF::Execution
