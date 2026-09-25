/// @file SF_globalDofResidual.cpp
/// @brief halo、GlobalDof 汇总与并行协调基础设施实现。

#include "SF_globalDofResidual.h"

#include "SF_MultiBlockMesh.h"
#include "backend/SF_communicationBackend.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace SF::Parallel {
namespace {

using MeshCommunication::HaloExchangePlan;
using MeshCommunication::HaloInterfacePoint;
using MeshCommunication::HaloInterfaceSyncGroup;

int blockOwner(const HaloExchangePlan& plan, int blockId, int size) {
    if (blockId < 0 || blockId >= (int)plan.blockInfos.size()) {
        throw std::runtime_error(
            "GlobalDof residual references an invalid block id.");
    }
    const int owner = blockId < (int)plan.blockOwnerRanks.size()
        ? plan.blockOwnerRanks[(size_t)blockId]
        : blockId % size;
    if (owner < 0 || owner >= size) {
        throw std::runtime_error(
            "GlobalDof residual references an invalid block owner rank.");
    }
    return owner;
}

void pointIJK(const Field& field, int interior, int& i, int& j, int& k) {
    const int pointCount = field.NX() * field.NY() * field.NZ();
    if (interior < 0 || interior >= pointCount) {
        throw std::runtime_error(
            "GlobalDof residual references an invalid interior point.");
    }
    i = interior % field.NX() + field.NG();
    j = (interior / field.NX()) % field.NY() + field.NG();
    k = interior / (field.NX() * field.NY()) + field.NG();
}

std::vector<int> displacements(const std::vector<int>& counts) {
    std::vector<int> result(counts.size(), 0);
    for (size_t rank = 1; rank < counts.size(); ++rank) {
        result[rank] = result[rank - 1] + counts[rank - 1];
    }
    return result;
}

int totalCount(const std::vector<int>& counts) {
    int result = 0;
    for (int count : counts) {
        if (count < 0) {
            throw std::runtime_error(
                "GlobalDof residual has a negative communication count.");
        }
        result += count;
    }
    return result;
}

std::vector<int> ownerNeighbours(
        const HaloExchangePlan& plan, int rank, int size) {
    std::set<int> neighbours;
    for (const HaloInterfaceSyncGroup& group : plan.interfaceSyncGroups) {
        if (group.ownerRank < 0 || group.ownerRank >= size) {
            throw std::runtime_error(
                "GlobalDof residual group has an invalid owner rank.");
        }
        for (const HaloInterfacePoint& point : group.points) {
            const int contributor = blockOwner(plan, point.blockId, size);
            if (rank == group.ownerRank && contributor != rank) {
                neighbours.insert(contributor);
            } else if (rank == contributor && group.ownerRank != rank) {
                neighbours.insert(group.ownerRank);
            }
        }
    }
    return {neighbours.begin(), neighbours.end()};
}

void exchange(
        const Backend::CommunicationBackend* backend,
        const std::vector<double>& send,
        const std::vector<int>& counts,
        const std::vector<int>& offsets,
        const std::vector<int>& neighbours,
        int rank,
        int tag,
        std::vector<double>& received) {
    if (counts[(size_t)rank] != (int)send.size()) {
        throw std::runtime_error(
            "GlobalDof residual payload differs from its static layout.");
    }
    std::copy(send.begin(), send.end(),
              received.begin() + offsets[(size_t)rank]);
    if (backend && backend->active()) {
        backend->exchangeNeighbours(
            send, counts, offsets, neighbours, tag, received);
    }
}

} // namespace

static void assembleImpl(
        const HaloExchangePlan& plan,
        std::vector<MeshBlockField>* blocks,
        const std::vector<Residual*>* residuals,
        Field* localField,
        Residual* localResidual,
        int localBlockId,
        const Backend::CommunicationBackend* backend) {
    if (plan.interfaceSyncGroups.empty()) return;
    if ((!blocks || blocks->empty()) && !localField) {
        throw std::runtime_error(
            "GlobalDof residual assembly requires local Field storage.");
    }

    const int rank = backend ? backend->rank() : 0;
    const int size = backend ? std::max(1, backend->size()) : 1;
    const int nVar = blocks && !blocks->empty()
        ? blocks->front().field.NVar() : localField->NVar();
    if (nVar <= 0) {
        throw std::runtime_error(
            "GlobalDof residual assembly requires conserved variables.");
    }
    if (blocks) {
        for (const MeshBlockField& block : *blocks) {
            if (block.field.NVar() != nVar) {
                throw std::runtime_error(
                    "GlobalDof residual blocks use different FluidStateModels.");
            }
        }
    }
    auto fieldForLocalBlock = [&](int blockId) -> Field& {
        if (blocks) {
            if (blockId < 0 || blockId >= (int)blocks->size()) {
                throw std::runtime_error(
                    "GlobalDof residual references unavailable block storage.");
            }
            return (*blocks)[(size_t)blockId].field;
        }
        if (!localField || blockId != localBlockId) {
            throw std::runtime_error(
                "GlobalDof residual requested a non-local Field.");
        }
        return *localField;
    };
    auto residualForLocalBlock = [&](int blockId) -> Residual& {
        if (residuals) {
            if (blockId < 0 || blockId >= (int)residuals->size()
                || !(*residuals)[(size_t)blockId]) {
                throw std::runtime_error("GlobalDof residual workspace is unavailable.");
            }
            return *(*residuals)[(size_t)blockId];
        }
        if (!localResidual || blockId != localBlockId) {
            throw std::runtime_error("GlobalDof residual requested a non-local workspace.");
        }
        return *localResidual;
    };

    const int contributionSize = nVar + 1;
    std::vector<int> contributionCounts((size_t)size, 0);
    for (const HaloInterfaceSyncGroup& group : plan.interfaceSyncGroups) {
        if (group.points.size() < 2) {
            throw std::runtime_error(
                "GlobalDof residual group must contain at least two replicas.");
        }
        for (const HaloInterfacePoint& point : group.points) {
            contributionCounts[(size_t)blockOwner(
                plan, point.blockId, size)] += contributionSize;
        }
    }
    const auto contributionOffsets = displacements(contributionCounts);
    std::vector<int> contributionCursor((size_t)size, 0);
    std::vector<int> pointOffsets;
    for (const auto& group : plan.interfaceSyncGroups) {
        pointOffsets.reserve(pointOffsets.size() + group.points.size());
        for (const auto& point : group.points) {
            const int owner = blockOwner(plan, point.blockId, size);
            pointOffsets.push_back(
                contributionOffsets[(size_t)owner]
                + contributionCursor[(size_t)owner]);
            contributionCursor[(size_t)owner] += contributionSize;
        }
    }

    std::vector<double> localContributions;
    localContributions.reserve((size_t)contributionCounts[(size_t)rank]);
    for (const auto& group : plan.interfaceSyncGroups) {
        for (const auto& point : group.points) {
            if (blockOwner(plan, point.blockId, size) != rank) continue;
            if (!std::isfinite(point.dualVolume)
                || point.dualVolume <= 1.0e-300) {
                throw std::runtime_error(
                    "GlobalDof residual has invalid dual-volume metadata.");
            }
            const Field& field = fieldForLocalBlock(point.blockId);
            const Residual& residual = residualForLocalBlock(point.blockId);
            int i = 0, j = 0, k = 0;
            pointIJK(field, point.interiorIndex, i, j, k);
            localContributions.push_back(point.dualVolume);
            for (int v = 0; v < nVar; ++v) {
                const double value = residual.local(i, j, k, v);
                if (!std::isfinite(value)) {
                    throw std::runtime_error(
                        "GlobalDof contributor produced a non-finite residual.");
                }
                localContributions.push_back(point.dualVolume * value);
            }
        }
    }

    const auto neighbours = ownerNeighbours(plan, rank, size);
    std::vector<double> allContributions(
        (size_t)totalCount(contributionCounts), 0.0);
    exchange(backend, localContributions, contributionCounts,
             contributionOffsets, neighbours, rank, 4111,
             allContributions);

    std::vector<std::vector<double>> ownerResults(
        plan.interfaceSyncGroups.size());
    size_t pointCursor = 0;
    for (size_t groupId = 0;
         groupId < plan.interfaceSyncGroups.size(); ++groupId) {
        const auto& group = plan.interfaceSyncGroups[groupId];
        if (group.ownerRank != rank) {
            pointCursor += group.points.size();
            continue;
        }
        double volume = 0.0;
        std::vector<double> integrated((size_t)nVar, 0.0);
        for (size_t point = 0; point < group.points.size(); ++point) {
            const size_t offset = (size_t)pointOffsets[pointCursor + point];
            volume += allContributions[offset];
            for (int v = 0; v < nVar; ++v) {
                integrated[(size_t)v] +=
                    allContributions[offset + (size_t)v + 1U];
            }
        }
        pointCursor += group.points.size();
        if (!std::isfinite(volume) || volume <= 1.0e-300) {
            throw std::runtime_error(
                "GlobalDof owner received a non-positive total dual volume.");
        }
        auto& result = ownerResults[groupId];
        result.resize((size_t)nVar);
        for (int v = 0; v < nVar; ++v) {
            result[(size_t)v] = integrated[(size_t)v] / volume;
            if (!std::isfinite(result[(size_t)v])) {
                throw std::runtime_error(
                    "GlobalDof owner assembled a non-finite residual.");
            }
        }
    }

    std::vector<int> resultCounts((size_t)size, 0);
    for (const auto& group : plan.interfaceSyncGroups) {
        resultCounts[(size_t)group.ownerRank] += nVar;
    }
    const auto resultOffsets = displacements(resultCounts);
    std::vector<int> resultCursor((size_t)size, 0);
    std::vector<int> groupOffsets(plan.interfaceSyncGroups.size(), -1);
    for (size_t groupId = 0;
         groupId < plan.interfaceSyncGroups.size(); ++groupId) {
        const int owner = plan.interfaceSyncGroups[groupId].ownerRank;
        groupOffsets[groupId] = resultOffsets[(size_t)owner]
            + resultCursor[(size_t)owner];
        resultCursor[(size_t)owner] += nVar;
    }

    std::vector<double> localResults;
    localResults.reserve((size_t)resultCounts[(size_t)rank]);
    for (size_t groupId = 0;
         groupId < plan.interfaceSyncGroups.size(); ++groupId) {
        if (plan.interfaceSyncGroups[groupId].ownerRank != rank) continue;
        localResults.insert(localResults.end(), ownerResults[groupId].begin(),
                            ownerResults[groupId].end());
    }
    std::vector<double> allResults(
        (size_t)totalCount(resultCounts), 0.0);
    exchange(backend, localResults, resultCounts, resultOffsets,
             neighbours, rank, 4112, allResults);

    for (size_t groupId = 0;
         groupId < plan.interfaceSyncGroups.size(); ++groupId) {
        const auto& group = plan.interfaceSyncGroups[groupId];
        const size_t offset = (size_t)groupOffsets[groupId];
        for (const auto& point : group.points) {
            if (blockOwner(plan, point.blockId, size) != rank) continue;
            Field& field = fieldForLocalBlock(point.blockId);
            Residual& residual = residualForLocalBlock(point.blockId);
            int i = 0, j = 0, k = 0;
            pointIJK(field, point.interiorIndex, i, j, k);
            for (int v = 0; v < nVar; ++v) {
                const double value = allResults[offset + (size_t)v];
                if (!std::isfinite(value)) {
                    throw std::runtime_error(
                        "GlobalDof residual broadcast is non-finite.");
                }
                residual.setGlobal(i, j, k, v, value);
            }
        }
    }
}

void GlobalDofResidualAssembler::assemble(
        const HaloExchangePlan& plan,
        std::vector<MeshBlockField>& blocks,
        const std::vector<Residual*>& residuals,
        const Backend::CommunicationBackend* backend) {
    assembleImpl(plan, &blocks, &residuals, nullptr, nullptr, -1, backend);
}

void GlobalDofResidualAssembler::assemble(
        const HaloExchangePlan& plan,
        Field& localField,
        Residual& residual,
        int localBlockId,
        const Backend::CommunicationBackend* backend) {
    assembleImpl(plan, nullptr, nullptr, &localField, &residual,
                 localBlockId, backend);
}

} // namespace SF::Parallel
