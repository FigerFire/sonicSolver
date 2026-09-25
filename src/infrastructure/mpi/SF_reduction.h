#pragma once

/// @file SF_reduction.h
/// @brief Reduction、Barrier 和 Broadcast 全局通信能力。

#include "backend/SF_communicationBackend.h"

namespace SF::Parallel {

class ReductionService {
public:
    explicit ReductionService(
        const Backend::CommunicationBackend* backend = nullptr)
        : backend_(backend) {}

    void attach(const Backend::CommunicationBackend* backend) {
        backend_ = backend;
    }
    double minimum(double value) const {
        return backend_ ? backend_->allReduceMin(value) : value;
    }
    double maximum(double value) const {
        return backend_ ? backend_->allReduceMax(value) : value;
    }
    double sum(double value) const {
        return backend_ ? backend_->allReduceSum(value) : value;
    }
    void sum(std::vector<double>& values) const {
        if (backend_) backend_->allReduceSum(values);
    }
    std::int64_t sum(std::int64_t value) const {
        return backend_ ? backend_->allReduceSum(value) : value;
    }
    bool allRanksAgree(bool value) const {
        return backend_ ? backend_->allRanksAgree(value) : value;
    }
    void broadcast(std::vector<double>& values, int root) const {
        if (backend_) backend_->broadcast(values, root);
    }
    void barrier() const {
        if (backend_) backend_->barrier();
    }

private:
    const Backend::CommunicationBackend* backend_ = nullptr;
};

} // namespace SF::Parallel
