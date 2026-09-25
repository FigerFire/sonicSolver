/// @file SF_executionRuntime.cpp
/// @brief Field freshness、halo 和装配契约的 ExecutionRuntime 实现。

#include "SF_executionRuntime.h"


#include <cmath>
#include <stdexcept>

namespace SF::Execution {

void Runtime::attachState(State::StateBundle& state) {
    state_ = &state;
    if (parallel_) parallel_->attachState(state);
}

void Runtime::prepare(const OperatorContract& contract) {
    contract.validate();
    ++statistics_.preparedOperators;
    if (!state_) {
        throw std::runtime_error(
            "ExecutionRuntime::prepare requires an attached StateBundle for '"
            + contract.name + "'.");
    }

    for (const auto& access : contract.accesses) {
        if (access.mode != AccessMode::ReadHalo) continue;
        ++statistics_.haloRequests;
        const bool needsExchange = state_->distributed.needsExchange(
            access.field, access.haloDepth);
        if (parallel_) {
            parallel_->ensureHalo(access.field, access.haloDepth);
        } else {
            state_->distributed.select(
                access.field, State::HaloSyncStage::None, access.haloDepth);
            state_->distributed.markSynchronized(
                access.field, access.haloDepth);
        }
        if (needsExchange) ++statistics_.haloExchanges;
        else ++statistics_.freshHaloReuses;
    }
}

void Runtime::finalize(const OperatorContract& contract) {
    contract.validate();
    ++statistics_.finalizedOperators;
    if (!state_) {
        throw std::runtime_error(
            "ExecutionRuntime::finalize requires an attached StateBundle for '"
            + contract.name + "'.");
    }

    bool assembleCanonicalFaces = false;
    std::vector<std::string> accumulatedDofs;
    std::vector<std::string> modified;
    for (const auto& access : contract.accesses) {
        if (access.mode == AccessMode::WriteOwned) {
            bool duplicate = false;
            for (const auto& name : modified) duplicate |= name == access.field;
            if (!duplicate) modified.push_back(access.field);
        } else if (access.mode == AccessMode::WriteCanonicalFace) {
            assembleCanonicalFaces = true;
        } else if (access.mode == AccessMode::AccumulateGlobalDof) {
            bool duplicate = false;
            for (const auto& name : accumulatedDofs) {
                duplicate |= name == access.field;
            }
            if (!duplicate) accumulatedDofs.push_back(access.field);
        }
    }

    for (const auto& name : modified) {
        if (parallel_) parallel_->markModified(name);
        else state_->distributed.markModified(name);
        ++statistics_.ownedWrites;
    }
    if (assembleCanonicalFaces || !accumulatedDofs.empty()) {
        throw std::runtime_error(
            "workspace synchronization must be invoked with explicit borrowed storage.");
    }
}

void Runtime::synchronizeCanonicalFaceFluxes(
        const std::vector<Field*>& fields,
        const std::vector<FluxField*>& fluxes,
        const std::vector<Residual*>& residuals) {
    if (!state_ || fields.size() != fluxes.size() || fields.size() != residuals.size()) {
        throw std::runtime_error("canonical flux workspace is not aligned with StateBundle patches.");
    }
    if (parallel_) parallel_->assembleCanonicalInterfaceFluxes(fields, fluxes, residuals);
    ++statistics_.canonicalFaceAssemblies;
}

void Runtime::accumulateGlobalDofResiduals(
        const std::vector<Field*>& fields,
        const std::vector<Residual*>& residuals) {
    if (!state_ || fields.size() != residuals.size()) {
        throw std::runtime_error("GlobalDof residual workspace is not aligned with StateBundle patches.");
    }
    if (parallel_) parallel_->accumulateGlobalDof("conservativeResidual", fields, residuals);
    ++statistics_.globalDofAssemblies;
}

double Runtime::globalMinimum(double value) {
    return parallel_
        ? parallel_->reduce(value, FDM::IParallelCoordinator::Reduction::Minimum)
        : value;
}

double Runtime::globalMaximum(double value) {
    return parallel_
        ? parallel_->reduce(value, FDM::IParallelCoordinator::Reduction::Maximum)
        : value;
}

double Runtime::globalSum(double value) {
    return parallel_
        ? parallel_->reduce(value, FDM::IParallelCoordinator::Reduction::Sum)
        : value;
}

void Runtime::globalSum(std::vector<double>& values) {
    if (parallel_) parallel_->reduceSum(values);
}

bool Runtime::ownsCanonicalEntity(std::int64_t entityId) const {
    if (entityId < 0) {
        throw std::runtime_error(
            "ExecutionRuntime received an invalid canonical entity id.");
    }
    if (!parallel_ || !parallel_->active()) return true;
    const std::int64_t ranks=static_cast<std::int64_t>(parallel_->size());
    if (ranks <= 0) {
        throw std::runtime_error(
            "ExecutionRuntime has an invalid parallel size.");
    }
    return entityId % ranks == static_cast<std::int64_t>(parallel_->rank());
}

bool Runtime::distributed() const {
    return parallel_ && parallel_->active();
}

void Runtime::copyCanonicalEntities(
        const std::vector<std::int64_t>& entityIds,
        std::vector<double>& values) {
    if (entityIds.size()!=values.size()) {
        throw std::runtime_error(
            "ExecutionRuntime canonical copy size mismatch.");
    }
    std::vector<unsigned char> ownerMask(entityIds.size(), 0);
    for (std::size_t n=0;n<entityIds.size();++n) {
        if (entityIds[n]<0 || !std::isfinite(values[n])) {
            throw std::runtime_error(
                "ExecutionRuntime canonical copy received invalid data.");
        }
        ownerMask[n] = ownsCanonicalEntity(entityIds[n]) ? 1 : 0;
        if (!ownerMask[n] && values[n]!=0.0) {
            throw std::runtime_error(
                "ExecutionRuntime canonical copy received a non-owner write.");
        }
    }
    if (parallel_) {
        parallel_->copyCanonicalEntities(entityIds, values, ownerMask);
    }
}

FDM::DistributedIndexRange Runtime::allocateDistributedIndices(
        std::int64_t localCount) {
    if (localCount < 0) {
        throw std::runtime_error(
            "ExecutionRuntime received a negative distributed size.");
    }
    return parallel_
        ? parallel_->distributeIndices(localCount)
        : FDM::DistributedIndexRange{0,localCount-1,localCount};
}

void Runtime::synchronizeTransient(
        const std::vector<State::DistributedFieldView>& fields) {
    if (parallel_) parallel_->synchronizeTransient(fields);
}

std::vector<int> Runtime::canonicalOwnerCells(
        Field& geometry,
        const std::vector<int>& candidates) {
    if (!parallel_) return candidates;
    std::vector<double> owner(
        static_cast<std::size_t>(geometry.TotalSize()),-1.0);
    for (int cell : candidates) {
        if (cell < 0 || cell >= geometry.TotalSize()) {
            throw std::runtime_error(
                "ExecutionRuntime owner selection received invalid cell.");
        }
        owner[static_cast<std::size_t>(cell)] =
            static_cast<double>(parallel_->rank());
    }
    auto view = State::workspaceView(
        "canonicalEulerianOwner",0,geometry,owner,1,geometry.NG(),
        State::HaloSyncStage::None,State::ExchangeKind::Identifier);
    parallel_->synchronizeTransient({view});
    std::vector<int> result;
    result.reserve(candidates.size());
    for (int cell : candidates) {
        const double remote = owner[static_cast<std::size_t>(cell)];
        if (!std::isfinite(remote) || remote < 0.0) {
            throw std::runtime_error(
                "ExecutionRuntime owner synchronization produced invalid rank.");
        }
        const int canonical = std::min(
            parallel_->rank(),static_cast<int>(std::llround(remote)));
        if (canonical == parallel_->rank()) result.push_back(cell);
    }
    return result;
}

} // namespace SF::Execution
