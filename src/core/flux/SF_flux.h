/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.07.01-----------*/

#pragma once

/// @file SF_flux.h
/// @brief 面通量存储容器。

#include <cstddef>
#include <vector>

namespace SF {

/// @brief cell-face拓扑上的守恒面通量。
///
/// FluxField以全局/分区局部faceId为索引保存`F_hat`。它是后续
/// owner-neighbour守恒装配的存储层，区别于Residual中按结构方向
/// 保存的旧半点通量缓存。
class FluxField {
public:
    /// @brief 默认构造,不分配内存。
    FluxField() = default;

    /// @brief 初始化面通量数组。
    /// @param faceCount 面数量。
    /// @param nVar 守恒变量个数。
    void setup(std::size_t faceCount, int nVar);

    /// @brief 清空所有面通量。
    void clear();

    /// @brief 面通量可写访问。
    /// @param faceId 面编号。
    /// @param var 守恒变量编号。
    double& operator()(std::size_t faceId, int var);

    /// @brief 面通量只读访问。
    double operator()(std::size_t faceId, int var) const;

    /// @brief 当前面数量。
    std::size_t faceCount() const { return faceCount_; }

    /// @brief 守恒变量数量。
    int variableCount() const { return nVar_; }

private:
    std::size_t faceCount_ = 0;
    int nVar_ = 0;
    std::vector<double> values_;

    std::size_t index(std::size_t faceId, int var) const;
};

} // namespace SF
