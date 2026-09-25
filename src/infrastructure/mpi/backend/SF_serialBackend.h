#pragma once

/// @file SF_serialBackend.h
/// @brief 单进程通信后端。

#include "SF_communicationBackend.h"

namespace SF::Parallel::Backend {

class SerialBackend final : public CommunicationBackend {
public:
    bool available() const override { return true; }
    bool active() const override { return false; }
    int rank() const override { return 0; }
    int size() const override { return 1; }

    bool allRanksAgree(bool localValue) const override { return localValue; }
    double allReduceMin(double localValue) const override { return localValue; }
    double allReduceMax(double localValue) const override { return localValue; }
    double allReduceSum(double localValue) const override { return localValue; }
    void allReduceSum(std::vector<double>&) const override {}
    void copyCanonical(
        const std::vector<std::int64_t>&,
        std::vector<double>&,
        const std::vector<unsigned char>&) const override {}
    std::int64_t exclusiveScanSum(std::int64_t) const override { return 0; }
    std::int64_t allReduceSum(std::int64_t localValue) const override {
        return localValue;
    }
    void broadcast(std::vector<double>&, int) const override {}
    void barrier() const override {}
    void exchangeNeighbours(
        const std::vector<double>&,
        const std::vector<int>&,
        const std::vector<int>&,
        const std::vector<int>&,
        int,
        std::vector<double>&) const override {}
};

} // namespace SF::Parallel::Backend
