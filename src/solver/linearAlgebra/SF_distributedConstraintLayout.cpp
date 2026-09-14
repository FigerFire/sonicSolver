/// @file SF_distributedConstraintLayout.cpp
/// @brief ConstraintGlobalDof owner row layout 的纯 infrastructure 实现。

#include "SF_distributedConstraintLayout.h"

#include <algorithm>
#include <stdexcept>
#include <string>

namespace SF::LinearAlgebra {

DistributedConstraintLayout DistributedConstraintLayout::build(
        const std::vector<GlobalConstraintDofId>& entities,
        FDM::IExecutionRuntime& runtime) {
    DistributedConstraintLayout result;
    if (entities.empty()) {
        throw std::runtime_error(
            "Distributed constraint layout cannot be built from an empty set.");
    }
    std::vector<std::int64_t> ids;
    ids.reserve(entities.size());
    for (const auto entity : entities) {
        if (!entity.valid()) {
            throw std::runtime_error(
                "Distributed constraint layout received an invalid entity.");
        }
        if (!result.rows_.emplace(entity.value(),-1).second) {
            throw std::runtime_error(
                "Distributed constraint layout received duplicate entity "
                +std::to_string(entity.value())+".");
        }
        ids.push_back(entity.value());
        if (runtime.ownsCanonicalEntity(entity.value())) {
            result.ownedEntities_.push_back(entity);
        }
    }

    // owner-local row 顺序必须由稳定 GlobalConstraintDof 决定，不能依赖
    // 各 partition 的局部遍历顺序；否则同一个约束在重分区后会获得不同 row。
    std::sort(
        result.ownedEntities_.begin(), result.ownedEntities_.end(),
        [](GlobalConstraintDofId a, GlobalConstraintDofId b) {
            return a.value() < b.value();
        });

    const FDM::DistributedIndexRange range=runtime.allocateDistributedIndices(
        static_cast<std::int64_t>(result.ownedEntities_.size()));
    result.firstRow_=range.first;
    result.lastRow_=range.last;
    result.globalSize_=range.total;
    if (result.globalSize_<=0) {
        throw std::runtime_error(
            "Distributed constraint layout has no global rows.");
    }

    std::vector<double> ownerRows(entities.size(),0.0);
    for (std::size_t n=0;n<entities.size();++n) {
        if (!runtime.ownsCanonicalEntity(entities[n].value())) continue;
        const auto found=std::lower_bound(
            result.ownedEntities_.begin(),result.ownedEntities_.end(),
            entities[n],[](GlobalConstraintDofId a,GlobalConstraintDofId b) {
                return a.value()<b.value();
            });
        if (found==result.ownedEntities_.end()
            || found->value()!=entities[n].value()) {
            throw std::runtime_error(
                "Distributed constraint layout lost an owned entity.");
        }
        const std::size_t owned=static_cast<std::size_t>(
            found-result.ownedEntities_.begin());
        ownerRows[n]=static_cast<double>(range.first
            +static_cast<std::int64_t>(owned)+1);
    }
    runtime.copyCanonicalEntities(ids,ownerRows);
    for (std::size_t n=0;n<entities.size();++n) {
        const double encoded=ownerRows[n];
        const auto rounded=static_cast<std::int64_t>(encoded);
        if (encoded<=0.0 || static_cast<double>(rounded)!=encoded
            || rounded-1<0 || rounded-1>=result.globalSize_) {
            throw std::runtime_error(
                "Distributed constraint layout received an invalid owner row.");
        }
        result.rows_[entities[n].value()]=rounded-1;
    }
    return result;
}

std::int64_t DistributedConstraintLayout::row(
        GlobalConstraintDofId entity) const {
    if (!entity.valid()) {
        throw std::runtime_error(
            "Distributed constraint row query received an invalid entity.");
    }
    const auto found=rows_.find(entity.value());
    if (found==rows_.end()) {
        throw std::runtime_error(
            "Distributed constraint layout does not contain entity "
            +std::to_string(entity.value())+".");
    }
    return found->second;
}

} // namespace SF::LinearAlgebra
