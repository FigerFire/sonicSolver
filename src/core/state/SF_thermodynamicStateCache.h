#pragma once

/// @file SF_thermodynamicStateCache.h
/// @brief EOS 派生热力学状态的记忆化缓存（Field 的 derived-data cache）。
///
/// 该对象只回答 "守恒状态 -> 热力学闭合结果如何缓存"，不持有守恒量本身。
/// Field 读取守恒量后调用 state()，缓存按 (cell, version) 失效。

#include "SF_fluidStateModel.h"

#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

namespace SF {
namespace State {

class ThermodynamicStateCache {
public:
    /// @brief 绑定提供热力学闭合的 FluidStateModel，并重置缓存与版本。
    void bind(const std::shared_ptr<const Physics::FluidStateModel::Model>& equations) {
        equations_ = equations;
        entries_.clear();
        version_ = 1;
    }

    /// @brief 设置哈希表容量（取 min(capacity, 4096)，与原逻辑一致）。
    void setCapacity(int capacity) {
        entries_.assign(
            static_cast<size_t>(std::min(capacity, 4096)), Entry{});
    }

    /// @brief 声明守恒时间层已写入，使缓存失效。
    void invalidate() {
        if (version_ == std::numeric_limits<std::uint64_t>::max()) {
            version_ = 1;
            for (auto& entry : entries_) {
                entry.cell = -1;
                entry.version = 0;
            }
            return;
        }
        ++version_;
    }

    /// @brief 清空缓存并回到初始版本。
    void reset() {
        entries_.clear();
        version_ = 1;
    }

    std::uint64_t version() const { return version_; }

    /// @brief 查询 (cell, version) 是否已有派生状态。
    /// @return 命中时返回缓存项指针；未命中返回 nullptr。
    const Physics::FluidStateModel::ThermodynamicState* find(
            int cell) const {
        if (entries_.empty()) return nullptr;
        const Entry& entry =
            entries_[static_cast<size_t>(cell) % entries_.size()];
        if (entry.cell != cell || entry.version != version_) return nullptr;
        return &entry.state;
    }

    /// @brief 计算未缓存的 EOS 闭合并写入缓存。
    ///
    /// 只有 cache miss 的调用者才会读取守恒状态，因此 Field 可以先查缓存，
    /// 避免命中路径仍然收集 q。
    Physics::FluidStateModel::ThermodynamicState computeAndStore(
            int cell, const double* q, int variableCount) const {
        requireClosure();
        auto state = equations_->close(q, variableCount);
        if (!entries_.empty()) {
            Entry& entry =
                entries_[static_cast<size_t>(cell) % entries_.size()];
            entry.cell = cell;
            entry.version = version_;
            entry.state = state;
        }
        return state;
    }

    /// @brief 以记忆化方式求取守恒状态 q 的 EOS 派生状态。
    Physics::FluidStateModel::ThermodynamicState state(
            const double* q, int variableCount, int cell) const {
        requireClosure();
        if (!entries_.empty()) {
            Entry& entry = entries_[
                static_cast<size_t>(cell) % entries_.size()];
            if (entry.cell == cell && entry.version == version_) {
                return entry.state;
            }
        }
        auto state = equations_->close(q, variableCount);
        if (!entries_.empty()) {
            Entry& entry = entries_[
                static_cast<size_t>(cell) % entries_.size()];
            entry.cell = cell;
            entry.version = version_;
            entry.state = state;
        }
        return state;
    }

private:
    void requireClosure() const {
        if (!equations_) {
            throw std::runtime_error(
                "ThermodynamicStateCache has no FluidStateModel closure.");
        }
    }

    struct Entry {
        int cell = -1;
        std::uint64_t version = 0;
        Physics::FluidStateModel::ThermodynamicState state;
    };

    mutable std::vector<Entry> entries_;
    std::uint64_t version_ = 1;
    std::shared_ptr<const Physics::FluidStateModel::Model> equations_;
};

} // namespace State
} // namespace SF
