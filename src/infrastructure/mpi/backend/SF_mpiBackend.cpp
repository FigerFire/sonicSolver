/// @file SF_mpiBackend.cpp
/// @brief MPI/HYPRE 的基础设施后端封装实现。

#include "SF_mpiBackend.h"

#include "core/interfaces/SF_log.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#if SF_USE_MPI
#include <mpi.h>
#endif

namespace SF::Parallel::Backend {

MPIBackend::MPIBackend(int& argc, char**& argv, bool requested)
    : requested_(requested) {
#if SF_USE_MPI
    available_ = true;
    if (!requested_) return;
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (!initialized) {
        MPI_Init(&argc, &argv);
        ownsMPI_ = true;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &size_);
    enabled_ = true;
    setLogOutputEnabled(rank_ == 0);
    SF::broadcast("MPI enabled: ", std::to_string(size_) + " rank(s)");
#else
    (void)argc;
    (void)argv;
    if (requested_) {
        SF::broadcast(
            "MPI warning: ", "MPI requested but this build has no MPI support.");
    }
#endif
}

MPIBackend::~MPIBackend() {
#if SF_USE_MPI
    if (ownsMPI_) {
        int finalized = 0;
        MPI_Finalized(&finalized);
        if (!finalized) MPI_Finalize();
    }
#endif
}

bool MPIBackend::allRanksAgree(bool localValue) const {
#if SF_USE_MPI
    if (!enabled_) return localValue;
    int input = localValue ? 1 : 0;
    int output = 0;
    MPI_Allreduce(&input, &output, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    return output != 0;
#else
    return localValue;
#endif
}

double MPIBackend::allReduceMin(double localValue) const {
#if SF_USE_MPI
    if (!enabled_) return localValue;
    double output = localValue;
    MPI_Allreduce(&localValue, &output, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    return output;
#else
    return localValue;
#endif
}

double MPIBackend::allReduceMax(double localValue) const {
#if SF_USE_MPI
    if (!enabled_) return localValue;
    double output = localValue;
    MPI_Allreduce(&localValue, &output, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return output;
#else
    return localValue;
#endif
}

double MPIBackend::allReduceSum(double localValue) const {
#if SF_USE_MPI
    if (!enabled_) return localValue;
    double output = 0.0;
    MPI_Allreduce(&localValue, &output, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return output;
#else
    return localValue;
#endif
}

void MPIBackend::allReduceSum(std::vector<double>& values) const {
#if SF_USE_MPI
    if (!enabled_ || values.empty()) return;
    std::vector<double> result(values.size(), 0.0);
    MPI_Allreduce(values.data(), result.data(), static_cast<int>(values.size()),
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    values.swap(result);
#else
    (void)values;
#endif
}

void MPIBackend::copyCanonical(
        const std::vector<std::int64_t>& entityIds,
        std::vector<double>& values,
        const std::vector<unsigned char>& ownerMask) const {
#if SF_USE_MPI
    if (!enabled_) return;
    const bool matchingLengths = entityIds.size() == values.size()
        && entityIds.size() == ownerMask.size();
    int globallyValid = matchingLengths ? 1 : 0;
    MPI_Allreduce(
        MPI_IN_PLACE, &globallyValid, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (globallyValid == 0) {
        throw std::runtime_error(
            "MPI canonical copy received mismatched id/value/owner lengths.");
    }
    std::unordered_set<std::int64_t> localIds;
    localIds.reserve(entityIds.size());
    std::vector<std::int64_t> sendIds;
    std::vector<double> sendValues;
    sendIds.reserve(entityIds.size() / static_cast<std::size_t>(size_) + 1);
    sendValues.reserve(sendIds.capacity());
    bool locallyValid = true;
    for (std::size_t n = 0; n < entityIds.size(); ++n) {
        const std::int64_t id = entityIds[n];
        const bool valid = id >= 0 && localIds.insert(id).second
            && std::isfinite(values[n]) && ownerMask[n] <= 1
            && (ownerMask[n] != 0 || values[n] == 0.0);
        locallyValid = locallyValid && valid;
    }
    globallyValid = locallyValid ? 1 : 0;
    MPI_Allreduce(
        MPI_IN_PLACE, &globallyValid, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (globallyValid == 0) {
        throw std::runtime_error(
            "MPI canonical copy received invalid or duplicate owner data.");
    }
    for (std::size_t n = 0; n < entityIds.size(); ++n) {
        if (ownerMask[n] != 0) {
            sendIds.push_back(entityIds[n]);
            sendValues.push_back(values[n]);
        }
    }
    const bool payloadFits = sendIds.size()
        <= static_cast<std::size_t>(std::numeric_limits<int>::max());
    globallyValid = payloadFits ? 1 : 0;
    MPI_Allreduce(
        MPI_IN_PLACE, &globallyValid, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (globallyValid == 0) {
        throw std::runtime_error(
            "MPI canonical copy payload exceeds MPI count range.");
    }
    const int localCount = static_cast<int>(sendIds.size());
    std::vector<int> counts(static_cast<std::size_t>(size_), 0);
    MPI_Allgather(&localCount, 1, MPI_INT, counts.data(), 1, MPI_INT,
                  MPI_COMM_WORLD);
    std::vector<int> displacements(counts.size(), 0);
    int total = 0;
    for (std::size_t rank = 0; rank < counts.size(); ++rank) {
        displacements[rank] = total;
        if (counts[rank] > std::numeric_limits<int>::max() - total)
            throw std::runtime_error(
                "MPI canonical copy payload exceeds MPI count range.");
        total += counts[rank];
    }
    std::vector<std::int64_t> allIds(static_cast<std::size_t>(total));
    std::vector<double> allValues(static_cast<std::size_t>(total));
    MPI_Allgatherv(sendIds.data(), localCount, MPI_INT64_T,
                   allIds.data(), counts.data(), displacements.data(),
                   MPI_INT64_T, MPI_COMM_WORLD);
    MPI_Allgatherv(sendValues.data(), localCount, MPI_DOUBLE,
                   allValues.data(), counts.data(), displacements.data(),
                   MPI_DOUBLE, MPI_COMM_WORLD);
    std::unordered_map<std::int64_t, double> received;
    received.reserve(allIds.size());
    for (std::size_t n = 0; n < allIds.size(); ++n) {
        const auto inserted = received.emplace(allIds[n], allValues[n]);
        if (!inserted.second) {
            throw std::runtime_error(
                "MPI canonical copy received duplicate owner entities.");
        }
    }
    bool complete = true;
    for (std::size_t n = 0; n < entityIds.size(); ++n) {
        const auto found = received.find(entityIds[n]);
        complete = complete && found != received.end()
            && std::isfinite(found->second);
    }
    globallyValid = complete ? 1 : 0;
    MPI_Allreduce(
        MPI_IN_PLACE, &globallyValid, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (globallyValid == 0) {
        throw std::runtime_error("MPI canonical copy is missing an owner value.");
    }
    for (std::size_t n = 0; n < entityIds.size(); ++n) {
        values[n] = received.at(entityIds[n]);
    }
#else
    (void)entityIds;
    (void)values;
    (void)ownerMask;
#endif
}

std::int64_t MPIBackend::exclusiveScanSum(std::int64_t localValue) const {
#if SF_USE_MPI
    if (!enabled_) return 0;
    std::int64_t offset = 0;
    MPI_Exscan(
        &localValue, &offset, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
    return rank_ == 0 ? 0 : offset;
#else
    (void)localValue;
    return 0;
#endif
}

std::int64_t MPIBackend::allReduceSum(std::int64_t localValue) const {
#if SF_USE_MPI
    if (!enabled_) return localValue;
    std::int64_t output = 0;
    MPI_Allreduce(
        &localValue, &output, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
    return output;
#else
    return localValue;
#endif
}

void MPIBackend::broadcast(std::vector<double>& values, int root) const {
#if SF_USE_MPI
    if (!enabled_) return;
    MPI_Bcast(
        values.data(), static_cast<int>(values.size()), MPI_DOUBLE,
        root, MPI_COMM_WORLD);
#else
    (void)values;
    (void)root;
#endif
}

void MPIBackend::barrier() const {
#if SF_USE_MPI
    if (enabled_) MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void MPIBackend::exchangeNeighbours(
        const std::vector<double>& send,
        const std::vector<int>& counts,
        const std::vector<int>& displacements,
        const std::vector<int>& neighbours,
        int tag,
        std::vector<double>& received) const {
#if SF_USE_MPI
    if (!enabled_) return;
    if (rank_ < 0 || rank_ >= static_cast<int>(counts.size())
        || counts.size() != displacements.size()
        || counts[static_cast<size_t>(rank_)] != static_cast<int>(send.size())) {
        throw std::runtime_error(
            "Neighbour payload layout differs from the local static plan.");
    }
    std::vector<MPI_Request> requests;
    requests.reserve(neighbours.size() * 2);
    for (int neighbour : neighbours) {
        requests.push_back(MPI_REQUEST_NULL);
        MPI_Irecv(
            received.data() + displacements[static_cast<size_t>(neighbour)],
            counts[static_cast<size_t>(neighbour)], MPI_DOUBLE, neighbour, tag,
            MPI_COMM_WORLD, &requests.back());
    }
    for (int neighbour : neighbours) {
        requests.push_back(MPI_REQUEST_NULL);
        MPI_Isend(
            send.data(), static_cast<int>(send.size()), MPI_DOUBLE,
            neighbour, tag, MPI_COMM_WORLD, &requests.back());
    }
    if (!requests.empty()) {
        MPI_Waitall(
            static_cast<int>(requests.size()), requests.data(),
            MPI_STATUSES_IGNORE);
    }
#else
    (void)send;
    (void)counts;
    (void)displacements;
    (void)neighbours;
    (void)tag;
    (void)received;
#endif
}

} // namespace SF::Parallel::Backend
