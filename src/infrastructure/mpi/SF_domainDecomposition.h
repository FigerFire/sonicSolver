#pragma once

/// @file SF_domainDecomposition.h
/// @brief 分布式自由度所有权和连续全局编号。

#include "core/interfaces/SF_executionRuntime.h"
#include "backend/SF_communicationBackend.h"

#include <stdexcept>

namespace SF::Parallel {

class DomainDecomposition {
public:
    explicit DomainDecomposition(
        const Backend::CommunicationBackend* backend = nullptr)
        : backend_(backend) {}

    void attach(const Backend::CommunicationBackend* backend) {
        backend_ = backend;
    }

    FDM::DistributedIndexRange distribute(std::int64_t localCount) const {
        if (localCount < 0) {
            throw std::runtime_error(
                "Domain decomposition received a negative local count.");
        }
        const std::int64_t first = backend_
            ? backend_->exclusiveScanSum(localCount) : 0;
        const std::int64_t total = backend_
            ? backend_->allReduceSum(localCount) : localCount;
        return {first, first + localCount - 1, total};
    }

private:
    const Backend::CommunicationBackend* backend_ = nullptr;
};

} // namespace SF::Parallel
