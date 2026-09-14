/// @file SF_globalDofSystem.cpp
/// @brief 以 GlobalDofId 装配守恒残差和稀疏线性行。

#include "SF_globalDofSystem.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace SF::LinearAlgebra {

GlobalDofId GlobalDofId::make(
        GlobalDofSpace space, std::int64_t entity, int component) {
    constexpr std::int64_t entityLimit = (std::int64_t{1} << 56);
    const auto spaceValue = static_cast<std::int64_t>(space);
    if (spaceValue <= 0 || spaceValue > 15
        || entity < 0 || entity >= entityLimit
        || component < 0 || component > 7) {
        throw std::runtime_error(
            "GlobalDofId received an invalid space/entity/component.");
    }
    return GlobalDofId(
        (spaceValue << 59)
        | (static_cast<std::int64_t>(component) << 56)
        | entity);
}

void GlobalDofRow::add(GlobalDofId column, double value) {
    if (!row_.valid() || !column.valid() || !std::isfinite(value)) {
        throw std::runtime_error(
            "GlobalDof row received an invalid identity/coefficient.");
    }
    for (auto& entry : coefficients_) {
        if (entry.column == column) {
            entry.value += value;
            if (!std::isfinite(entry.value)) {
                throw std::runtime_error(
                    "GlobalDof row coefficient accumulation is non-finite.");
            }
            return;
        }
    }
    coefficients_.push_back({column, value});
}

StaticDistributedNumbering::StaticDistributedNumbering(
        std::int64_t firstRow,
        std::int64_t lastRow,
        std::int64_t globalSize,
        const std::vector<std::pair<GlobalDofId, std::int64_t>>& entries)
    : firstRow_(firstRow), lastRow_(lastRow), globalSize_(globalSize) {
    if (globalSize_ <= 0 || firstRow_ < 0 || lastRow_ < firstRow_
        || lastRow_ >= globalSize_) {
        throw std::runtime_error(
            "Distributed numbering has invalid row bounds.");
    }
    for (const auto& entry : entries) {
        if (!entry.first.valid() || entry.second < 0
            || entry.second >= globalSize_) {
            throw std::runtime_error(
                "Distributed numbering contains an invalid mapping.");
        }
        const auto inserted = rows_.emplace(
            entry.first.value(), entry.second);
        if (!inserted.second && inserted.first->second != entry.second) {
            throw std::runtime_error(
                "One GlobalDof maps to multiple distributed rows.");
        }
    }
}

std::int64_t StaticDistributedNumbering::row(GlobalDofId dof) const {
    const auto found = rows_.find(dof.value());
    if (!dof.valid() || found == rows_.end()) {
        throw std::runtime_error(
            "Distributed numbering does not contain GlobalDof "
            + std::to_string(dof.value()) + ".");
    }
    return found->second;
}

DistributedLinearSystem::DistributedLinearSystem(
        FDM::LinearSolverConfig config,
        const IDistributedNumbering& numbering)
    : numbering_(numbering), session_(std::move(config)) {}

SolveResult DistributedLinearSystem::solve(
        const GlobalDofSystem& equations) {
    const std::int64_t localSize =
        numbering_.lastRow() - numbering_.firstRow() + 1;
    if ((std::int64_t)equations.rows.size() != localSize) {
        throw std::runtime_error(
            "GlobalDof equation count does not match owned row range.");
    }
    SparseSystem backend;
    backend.firstRow = numbering_.firstRow();
    backend.lastRow = numbering_.lastRow();
    backend.globalSize = numbering_.globalSize();
    backend.rows.reserve(equations.rows.size());
    backend.rhs.reserve(equations.rows.size());
    backend.initialGuess.reserve(equations.rows.size());
    for (size_t local = 0; local < equations.rows.size(); ++local) {
        const GlobalDofRow& source = equations.rows[local];
        SparseRow target;
        target.globalRow = numbering_.row(source.row());
        const std::int64_t expected = backend.firstRow
            + (std::int64_t)local;
        if (target.globalRow != expected) {
            throw std::runtime_error(
                "GlobalDof rows are not ordered by the owned distributed "
                "row interval.");
        }
        for (const auto& coefficient : source.coefficients()) {
            target.columns.push_back(numbering_.row(coefficient.column));
            target.values.push_back(coefficient.value);
        }
        backend.rows.push_back(std::move(target));
        backend.rhs.push_back(source.rightHandSide());
        backend.initialGuess.push_back(source.initialGuess());
        backend.blockRoles.push_back(static_cast<int>(source.row().space()));
    }
    return session_.solve(backend);
}

void DistributedLinearSystem::invalidateStructure() {
    session_.invalidateStructure();
}

const ReuseStatistics& DistributedLinearSystem::statistics() const {
    return session_.statistics();
}

} // namespace SF::LinearAlgebra
