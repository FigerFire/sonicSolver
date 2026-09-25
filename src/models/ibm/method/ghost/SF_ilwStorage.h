/*--------------Sonic Fluid-------------------*/

#pragma once

/// @file SF_ilwStorage.h
/// @brief IBM-ILW 专属几何样本和预计算闭合计划存储。
///
/// 这些数据是 IBM ghost-cell ILW 的辅助几何状态，而不是通用守恒 Field
/// 状态。存储由 IBM 前处理对象拥有，不进入通用 Field。

#include "SF_ilwPlan.h"

#include <algorithm>
#include <cstddef>
#include <vector>

namespace SF::IBM::GhostILW {

/// @brief 每个 IBM patch 的 ILW 样本索引和预计算拟合计划。
class Storage {
public:
    void setup(std::size_t cellCount) {
        fluidOffset_.assign(cellCount, 0); fluidCount_.assign(cellCount, 0);
        normalOffset_.assign(cellCount, 0); normalCount_.assign(cellCount, 0);
        planIndex_.assign(cellCount, -1);
        fluidSamples_.clear(); normalSamples_.clear(); plans_.clear();
    }
    void clearSamples() {
        std::fill(fluidOffset_.begin(), fluidOffset_.end(), 0);
        std::fill(fluidCount_.begin(), fluidCount_.end(), 0);
        std::fill(normalOffset_.begin(), normalOffset_.end(), 0);
        std::fill(normalCount_.begin(), normalCount_.end(), 0);
        fluidSamples_.clear(); normalSamples_.clear();
    }
    void clearPlans() { std::fill(planIndex_.begin(), planIndex_.end(), -1); plans_.clear(); }
    void setFluidSamples(std::size_t cell, const std::vector<int>& values) {
        fluidOffset_[cell] = static_cast<int>(fluidSamples_.size());
        fluidCount_[cell] = static_cast<int>(values.size());
        fluidSamples_.insert(fluidSamples_.end(), values.begin(), values.end());
    }
    void setNormalSamples(std::size_t cell, const std::vector<int>& values) {
        normalOffset_[cell] = static_cast<int>(normalSamples_.size());
        normalCount_[cell] = static_cast<int>(values.size());
        normalSamples_.insert(normalSamples_.end(), values.begin(), values.end());
    }
    int fluidCount(std::size_t cell) const { return fluidCount_[cell]; }
    int fluidSample(std::size_t cell, int n) const { return fluidSamples_[static_cast<std::size_t>(fluidOffset_[cell] + n)]; }
    int normalCount(std::size_t cell) const { return normalCount_[cell]; }
    int normalSample(std::size_t cell, int n) const { return normalSamples_[static_cast<std::size_t>(normalOffset_[cell] + n)]; }
    void setPlan(std::size_t cell, const IBMILWPointPlan& plan) {
        planIndex_[cell] = static_cast<int>(plans_.size()); plans_.push_back(plan);
    }
    bool hasPlan(std::size_t cell) const {
        const int index = planIndex_[cell];
        return index >= 0 && index < static_cast<int>(plans_.size()) && plans_[static_cast<std::size_t>(index)].valid;
    }
    const IBMILWPointPlan& plan(std::size_t cell) const {
        static const IBMILWPointPlan invalidPlan;
        const int index = planIndex_[cell];
        return (index >= 0 && index < static_cast<int>(plans_.size())) ? plans_[static_cast<std::size_t>(index)] : invalidPlan;
    }
private:
    std::vector<int> fluidOffset_, fluidCount_, normalOffset_, normalCount_;
    std::vector<int> fluidSamples_, normalSamples_, planIndex_;
    std::vector<IBMILWPointPlan> plans_;
};

} // namespace SF::IBM::GhostILW
