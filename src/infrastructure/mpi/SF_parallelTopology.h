#pragma once

/// @file SF_parallelTopology.h
/// @brief 从 mesh communication plan 派生的执行单元拓扑。

#include "SF_communicationPlan.h"

#include <algorithm>
#include <vector>

namespace SF::Parallel {

class ParallelTopology {
public:
    void configure(
        const MeshCommunication::HaloExchangePlan* plan,
        int localBlockId,
        int rank,
        int size) {
        plan_ = plan;
        localBlockId_ = localBlockId;
        rank_ = rank;
        size_ = std::max(1, size);
        processorInterfaces_.clear();
        if (!plan_) return;
        for (const auto& interface : plan_->interfaces) {
            const int owner = ownerRank(interface.ownerBlock);
            const int neighbour = ownerRank(interface.neighbourBlock);
            if (owner == neighbour) continue;
            MeshCommunication::ProcessorInterface processor;
            static_cast<MeshCommunication::MeshInterface&>(processor) =
                interface;
            processor.ownerRank = owner;
            processor.neighbourRank = neighbour;
            processorInterfaces_.push_back(processor);
        }
    }

    const MeshCommunication::HaloExchangePlan* plan() const {
        return plan_;
    }
    int localBlockId() const { return localBlockId_; }
    int rank() const { return rank_; }
    int size() const { return size_; }
    const std::vector<MeshCommunication::ProcessorInterface>&
    processorInterfaces() const {
        return processorInterfaces_;
    }

    int ownerRank(int blockId) const {
        if (blockId < 0) return -1;
        if (plan_ && blockId < static_cast<int>(plan_->blockOwnerRanks.size())) {
            return plan_->blockOwnerRanks[static_cast<size_t>(blockId)];
        }
        return blockId % size_;
    }

private:
    const MeshCommunication::HaloExchangePlan* plan_ = nullptr;
    int localBlockId_ = -1;
    int rank_ = 0;
    int size_ = 1;
    std::vector<MeshCommunication::ProcessorInterface> processorInterfaces_;
};

} // namespace SF::Parallel
