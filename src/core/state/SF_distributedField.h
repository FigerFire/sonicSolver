#pragma once

/// @file SF_distributedField.h
/// @brief 分布式状态的非拥有字段视图和统一交换语义。

#include "core/field/SF_field.h"
#include "SF_scalarField.h"
#include "SF_stateLayout.h"

#include <cstddef>
#include <functional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace SF::State {

/// @brief 字段在 patch 接口上的交换语义。
enum class ExchangeKind {
    State,              ///< 普通状态：填充 halo，并统一共享网格点。
    Identifier,         ///< 离散编号：只复制 donor，不进行数值平均。
    CanonicalFaceFlux,  ///< owner 产生的唯一面通量。
    None                ///< 本地派生量，不参与分布式交换。
};

/// @brief 任意 cell/face 状态的非拥有、可读写分布式视图。
///
/// 数据的所有权始终留在 Field、ScalarField 或算法 workspace 中。并行层只消费
/// 此视图，不需要知道字段属于 Level Set、湍流还是某个多相模型。
struct DistributedFieldView {
    std::string name;
    int blockId = -1;
    Field* geometry = nullptr;
    FieldLocation location = FieldLocation::Cell;
    ExchangeKind exchange = ExchangeKind::State;
    int components = 1;
    // 过渡字段：表示存储可提供的最大 halo capacity，不是某个格式的需求。
    // 新算子的实际深度只写在 Execution::OperatorContract::ReadHalo 中。
    int haloDepth = 0;
    HaloSyncStage stages = HaloSyncStage::None;
    std::function<double(int, int)> read;
    std::function<void(int, int, double)> write;

    /// @brief 检查视图尺寸和访问器是否完整。
    void validate() const {
        if (name.empty() || !geometry || components <= 0
            || haloDepth < 0 || !read || !write) {
            throw std::runtime_error(
                "DistributedFieldView '" + name
                + "' has an incomplete storage contract.");
        }
        if (haloDepth > geometry->NG()) {
            throw std::runtime_error(
                "DistributedFieldView '" + name
                + "' requests halo depth larger than the allocated Field halo.");
        }
    }
};

/// @brief 为主守恒 Field 创建多分量视图。
inline DistributedFieldView conservativeView(
        std::string name, int blockId, Field& field,
        int haloDepth, HaloSyncStage stages) {
    DistributedFieldView view;
    view.name = std::move(name);
    view.blockId = blockId;
    view.geometry = &field;
    view.components = field.NVar();
    view.haloDepth = haloDepth;
    view.stages = stages;
    view.read = [&field](int cell, int component) {
        int i = 0, j = 0, k = 0;
        field.getIJK(cell, i, j, k);
        return field(i, j, k, component);
    };
    view.write = [&field](int cell, int component, double value) {
        int i = 0, j = 0, k = 0;
        field.getIJK(cell, i, j, k);
        field(i, j, k, component) = value;
    };
    return view;
}

/// @brief 注册主守恒存储，不在注册点写入任何格式相关 halo 深度或 stage。
inline DistributedFieldView conservativeView(
        std::string name, int blockId, Field& field) {
    return conservativeView(
        std::move(name), blockId, field,
        field.NG(), HaloSyncStage::None);
}

/// @brief 为 ScalarField 创建单分量视图。
inline DistributedFieldView scalarView(
        std::string name, int blockId, Field& geometry,
        ScalarField& field, int haloDepth, HaloSyncStage stages,
        ExchangeKind exchange = ExchangeKind::State) {
    if (!field.isCompatibleWith(geometry)) {
        throw std::runtime_error(
            "Distributed scalar '" + name + "' does not match its geometry.");
    }
    DistributedFieldView view;
    view.name = std::move(name);
    view.blockId = blockId;
    view.geometry = &geometry;
    view.exchange = exchange;
    view.haloDepth = haloDepth;
    view.stages = stages;
    view.read = [&field](int cell, int) { return field.values().at((size_t)cell); };
    view.write = [&field](int cell, int, double value) {
        field.values().at((size_t)cell) = value;
    };
    return view;
}

/// @brief 注册标量存储；halo capacity 来自网格，实际需求由 operator contract 给出。
inline DistributedFieldView scalarView(
        std::string name, int blockId, Field& geometry,
        ScalarField& field,
        ExchangeKind exchange = ExchangeKind::State) {
    return scalarView(
        std::move(name), blockId, geometry, field,
        geometry.NG(), HaloSyncStage::None, exchange);
}

/// @brief 为按 component-major 存储的算法数组创建视图。
inline DistributedFieldView workspaceView(
        std::string name, int blockId, Field& geometry,
        std::vector<double>& values, int components,
        int haloDepth, HaloSyncStage stages,
        ExchangeKind exchange = ExchangeKind::State,
        FieldLocation location = FieldLocation::Cell) {
    const size_t expected = (size_t)geometry.TotalSize() * (size_t)components;
    if (components <= 0 || values.size() != expected) {
        throw std::runtime_error(
            "Distributed workspace '" + name + "' has an invalid size.");
    }
    DistributedFieldView view;
    view.name = std::move(name);
    view.blockId = blockId;
    view.geometry = &geometry;
    view.location = location;
    view.exchange = exchange;
    view.components = components;
    view.haloDepth = haloDepth;
    view.stages = stages;
    const int total = geometry.TotalSize();
    view.read = [&values, total](int cell, int component) {
        return values.at((size_t)component * (size_t)total + (size_t)cell);
    };
    view.write = [&values, total](int cell, int component, double value) {
        values.at((size_t)component * (size_t)total + (size_t)cell) = value;
    };
    return view;
}

/// @brief 注册算法工作区存储；只记录可用容量，不声明任何格式或同步 stage。
inline DistributedFieldView workspaceView(
        std::string name, int blockId, Field& geometry,
        std::vector<double>& values, int components,
        ExchangeKind exchange = ExchangeKind::State,
        FieldLocation location = FieldLocation::Cell) {
    return workspaceView(
        std::move(name), blockId, geometry, values, components,
        geometry.NG(), HaloSyncStage::None, exchange, location);
}

/// @brief 把若干独立 ScalarField 组合为一个多分量 cell field。
inline DistributedFieldView scalarComponentsView(
        std::string name, int blockId, Field& geometry,
        std::vector<ScalarField*> components,
        int haloDepth, HaloSyncStage stages,
        ExchangeKind exchange = ExchangeKind::State) {
    if (components.empty()) {
        throw std::runtime_error(
            "Distributed component field '" + name + "' is empty.");
    }
    for (const auto* component : components) {
        if (!component || !component->isCompatibleWith(geometry)) {
            throw std::runtime_error(
                "Distributed component field '" + name
                + "' does not match its geometry.");
        }
    }
    DistributedFieldView view;
    view.name = std::move(name);
    view.blockId = blockId;
    view.geometry = &geometry;
    view.exchange = exchange;
    view.components = (int)components.size();
    view.haloDepth = haloDepth;
    view.stages = stages;
    view.read = [components](int cell, int component) {
        return components.at((size_t)component)->values().at((size_t)cell);
    };
    view.write = [components](int cell, int component, double value) {
        components.at((size_t)component)->values().at((size_t)cell) = value;
    };
    return view;
}

/// @brief 一个状态包拥有的全部分布式字段注册表。
class DistributedFieldRegistry {
public:
    void add(DistributedFieldView view) {
        view.validate();
        for (const auto& item : fields_) {
            if (item.view.name == view.name
                && item.view.blockId == view.blockId) {
                throw std::runtime_error(
                    "Distributed field '" + view.name
                    + "' is already registered for block "
                    + std::to_string(view.blockId) + ".");
            }
        }
        fields_.push_back({std::move(view), 1, 0, 0});
    }

    void clear() { fields_.clear(); }
    bool empty() const { return fields_.empty(); }
    size_t size() const { return fields_.size(); }

    void markModified(const std::string& name) {
        bool found = false;
        for (auto& item : fields_) {
            if (item.view.name == name) {
                ++item.generation;
                item.availableDepth = 0;
                found = true;
            }
        }
        if (!found) {
            throw std::runtime_error(
                "Cannot mark unregistered distributed field '" + name + "'.");
        }
    }

    std::vector<DistributedFieldView*> select(
            const std::string& name, HaloSyncStage stage,
            int requiredDepth = 0) {
        std::vector<DistributedFieldView*> selected;
        for (auto& item : fields_) {
            if (item.view.name != name) continue;
            if (stage != HaloSyncStage::None
                && !contains(item.view.stages, stage)) continue;
            if (requiredDepth > item.view.haloDepth) {
                throw std::runtime_error(
                    "Distributed field '" + name
                    + "' cannot provide the requested halo depth "
                    + std::to_string(requiredDepth) + ".");
            }
            selected.push_back(&item.view);
        }
        return selected;
    }

    std::vector<std::string> names(HaloSyncStage stage) const {
        std::vector<std::string> result;
        for (const auto& item : fields_) {
            if (!contains(item.view.stages, stage)) continue;
            bool duplicate = false;
            for (const auto& name : result) duplicate |= name == item.view.name;
            if (!duplicate) result.push_back(item.view.name);
        }
        return result;
    }

    bool needsExchange(const std::string& name, int depth) const {
        bool found = false;
        for (const auto& item : fields_) {
            if (item.view.name != name) continue;
            found = true;
            if (item.synchronizedGeneration != item.generation
                || item.availableDepth < depth) return true;
        }
        if (!found) {
            throw std::runtime_error(
                "Requested halo for unregistered distributed field '"
                + name + "'.");
        }
        return false;
    }

    void markSynchronized(const std::string& name, int depth) {
        for (auto& item : fields_) {
            if (item.view.name != name) continue;
            item.synchronizedGeneration = item.generation;
            if (depth > item.availableDepth) item.availableDepth = depth;
        }
    }

private:
    struct Entry {
        DistributedFieldView view;
        unsigned long long generation = 1;
        unsigned long long synchronizedGeneration = 0;
        int availableDepth = 0;
    };
    std::vector<Entry> fields_;
};

} // namespace SF::State
