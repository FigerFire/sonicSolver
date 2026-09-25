#pragma once

/// @file SF_mpiBackend.h
/// @brief MPI 通信后端；头文件不暴露 mpi.h。

#include "SF_communicationBackend.h"

namespace SF::Parallel::Backend {

class MPIBackend final : public CommunicationBackend {
public:
    MPIBackend(int& argc, char**& argv, bool requested);
    ~MPIBackend() override;

    MPIBackend(const MPIBackend&) = delete;
    MPIBackend& operator=(const MPIBackend&) = delete;

    bool available() const override { return available_; }
    bool active() const override { return enabled_; }
    int rank() const override { return rank_; }
    int size() const override { return size_; }

    bool allRanksAgree(bool localValue) const override;
    double allReduceMin(double localValue) const override;
    double allReduceMax(double localValue) const override;
    double allReduceSum(double localValue) const override;
    void allReduceSum(std::vector<double>& values) const override;
    void copyCanonical(
        const std::vector<std::int64_t>& entityIds,
        std::vector<double>& values,
        const std::vector<unsigned char>& ownerMask) const override;
    std::int64_t exclusiveScanSum(std::int64_t localValue) const override;
    std::int64_t allReduceSum(std::int64_t localValue) const override;
    void broadcast(std::vector<double>& values, int root) const override;
    void barrier() const override;
    void exchangeNeighbours(
        const std::vector<double>& send,
        const std::vector<int>& counts,
        const std::vector<int>& displacements,
        const std::vector<int>& neighbours,
        int tag,
        std::vector<double>& received) const override;

private:
    bool requested_ = false;
    bool available_ = false;
    bool enabled_ = false;
    bool ownsMPI_ = false;
    int rank_ = 0;
    int size_ = 1;
};

} // namespace SF::Parallel::Backend
