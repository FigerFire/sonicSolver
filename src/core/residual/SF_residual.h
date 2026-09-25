/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.07.01-----------*/

#pragma once

/// @file SF_residual.h
/// @brief 守恒残差与源项存储容器。

#include <cstddef>
#include <vector>

namespace SF {

/// @brief 结构网格守恒残差容器。
///
/// Residual只保存数值通量差更新所需的面通量缓存和显式源项:
/// `resX/resY/resZ` 对应计算坐标三个方向的半点面通量，
/// `source` 对应单元/点上的显式源项。索引布局与Field守恒量一致:
/// `[k][j][i][v]`。
class Residual {
public:
    /// @brief 默认构造,不分配内存。
    Residual() = default;

    /// @brief 初始化残差存储。
    /// @param mx 含虚胞的x方向存储点数。
    /// @param my 含虚胞的y方向存储点数。
    /// @param mz 含虚胞的z方向存储点数。
    /// @param nVar 守恒变量个数。
    void setup(int mx, int my, int mz, int nVar);

    /// @brief 查询既有 allocation 是否仍可服务指定 Field 形状；不修改存储。
    bool isCompatibleWith(int mx, int my, int mz, int nVar) const {
        return mx_ == mx && my_ == my && mz_ == mz && nVar_ == nVar;
    }

    /// @brief 清空面通量残差和源项。
    void clear();

    /// @brief 只清空源项。
    void clearSource();

    /// @brief x方向半点面通量缓存。
    double& x(int i, int j, int k, int v);
    double x(int i, int j, int k, int v) const;

    /// @brief y方向半点面通量缓存。
    double& y(int i, int j, int k, int v);
    double y(int i, int j, int k, int v) const;

    /// @brief z方向半点面通量缓存。
    double& z(int i, int j, int k, int v);
    double z(int i, int j, int k, int v) const;

    /// @brief 显式源项缓存。
    double& source(int i, int j, int k, int v);
    double source(int i, int j, int k, int v) const;

    /// @brief x方向残差数组长度,用于兼容旧诊断。
    std::size_t xSize() const { return resX_.size(); }

    /// @brief 保存当前算子完整装配后的本地强形式残差。
    void stageLocal(int i, int j, int k, int v, double value);
    /// @brief 读取本地强形式残差，供 GlobalDof 积分装配使用。
    double local(int i, int j, int k, int v) const;
    /// @brief 写入按对偶体积汇总后的 GlobalDof 强形式残差。
    void setGlobal(int i, int j, int k, int v, double value);
    /// @brief 当前点是否已有 Runtime 汇总后的 GlobalDof 残差。
    bool hasGlobal(int i, int j, int k) const;
    /// @brief 读取 Runtime 汇总后的 GlobalDof 残差。
    double global(int i, int j, int k, int v) const;
    /// @brief 清除上一算子阶段的 GlobalDof 残差标记。
    void clearGlobal();

private:
    int mx_ = 0;
    int my_ = 0;
    int mz_ = 0;
    int nVar_ = 0;
    std::vector<double> resX_;
    std::vector<double> resY_;
    std::vector<double> resZ_;
    std::vector<double> source_;
    std::vector<double> local_;
    std::vector<double> global_;
    std::vector<unsigned char> globalMask_;

    /// @brief 计算平铺数组索引。
    std::size_t index(int i, int j, int k, int v) const;
    std::size_t pointIndex(int i, int j, int k) const;
};

} // namespace SF
