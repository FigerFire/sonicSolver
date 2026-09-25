/// @file SF_localExecutionRuntime.cpp
/// @brief Serial execution-contract implementation.

#include "SF_localExecutionRuntime.h"

namespace SF::Execution {

void LocalRuntime::prepare(const OperatorContract& contract) {
    contract.validate();
    if (!state_) return;
    for (const auto& access : contract.accesses) {
        if (access.mode != AccessMode::ReadHalo) continue;
        state_->distributed.select(
            access.field,State::HaloSyncStage::None,access.haloDepth);
        state_->distributed.markSynchronized(access.field,access.haloDepth);
    }
}

void LocalRuntime::finalize(const OperatorContract& contract) {
    contract.validate();
    if (!state_) return;
    std::vector<std::string> modified;
    for (const auto& access : contract.accesses) {
        if (access.mode != AccessMode::WriteOwned) continue;
        bool duplicate = false;
        for (const auto& name : modified) duplicate |= name == access.field;
        if (duplicate) continue;
        state_->distributed.markModified(access.field);
        modified.push_back(access.field);
    }
}

bool LocalRuntime::ownsCanonicalEntity(std::int64_t entityId) const {
    if (entityId < 0) {
        throw std::runtime_error("Local runtime received an invalid canonical entity id.");
    }
    return true;
}

void LocalRuntime::copyCanonicalEntities(
        const std::vector<std::int64_t>& entityIds,
        std::vector<double>& values) {
    if (entityIds.size() != values.size()) {
        throw std::runtime_error("Local runtime canonical copy size mismatch.");
    }
    for (std::size_t n = 0; n < entityIds.size(); ++n) {
        if (entityIds[n] < 0 || !std::isfinite(values[n])) {
            throw std::runtime_error("Local runtime received invalid canonical copy data.");
        }
    }
}

FDM::DistributedIndexRange LocalRuntime::allocateDistributedIndices(
        std::int64_t localCount) {
    if (localCount < 0) {
        throw std::runtime_error("Local runtime received a negative distributed size.");
    }
    return {0,localCount-1,localCount};
}

} // namespace SF::Execution
