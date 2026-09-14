/// @file SF_distributedRowMap.cpp
/// @brief GlobalDofId 到 owner 与分布式线性行号的独立映射。

#include "SF_distributedRowMap.h"

#include "SF_distributedField.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace SF::LinearAlgebra {

DistributedRowMap makeDistributedRowMap(
        const Field& field,
        FDM::IExecutionRuntime* runtime,
        const CellSelector& selected) {
    if (!selected) {
        throw std::runtime_error(
            "Distributed row-map requires an explicit cell selector.");
    }

    DistributedRowMap map;
    map.global.assign((size_t)field.TotalSize(), -1);
    const int ng = field.NG();
    std::vector<int> candidates;
    for (int k = ng; k < ng + field.NZ(); ++k) {
        for (int j = ng; j < ng + field.NY(); ++j) {
            for (int i = ng; i < ng + field.NX(); ++i) {
                if (selected(i, j, k)) {
                    candidates.push_back(field.getIdx(i, j, k));
                }
            }
        }
    }

    map.localCells = runtime
        ? runtime->canonicalOwnerCells(
            const_cast<Field&>(field),candidates)
        : std::move(candidates);

    const std::int64_t local = (std::int64_t)map.localCells.size();
    const FDM::DistributedIndexRange ownership = runtime
        ? runtime->allocateDistributedIndices(local)
        : FDM::DistributedIndexRange{0, local - 1, local};
    map.first = ownership.first;
    map.last = ownership.first + local - 1;
    map.total = ownership.total;
    for (size_t n = 0; n < map.localCells.size(); ++n) {
        map.global[(size_t)map.localCells[n]] =
            ownership.first + (std::int64_t)n;
    }

    std::vector<double> encoded((size_t)field.TotalSize(), -1.0);
    for (size_t n = 0; n < map.localCells.size(); ++n) {
        encoded[(size_t)map.localCells[n]] =
            (double)(ownership.first + (std::int64_t)n);
    }
    if (runtime) {
        auto view = State::workspaceView(
            "legacyDistributedGlobalRow", 0, const_cast<Field&>(field),
            encoded, 1, field.NG(), State::HaloSyncStage::None,
            State::ExchangeKind::Identifier);
        runtime->synchronizeTransient({view});
    }
    for (size_t n = 0; n < encoded.size(); ++n) {
        if (encoded[n] < 0.0) continue;
        const double rounded = std::round(encoded[n]);
        if (std::abs(encoded[n] - rounded) > 1.0e-9
            || rounded > 9007199254740992.0) {
            throw std::runtime_error(
                "Distributed row-map halo id is invalid.");
        }
        if (map.global[n] < 0) map.global[n] = (std::int64_t)rounded;
    }
    return map;
}

} // namespace SF::LinearAlgebra
