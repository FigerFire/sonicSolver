/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.07.02-----------*/

#pragma once

/// @file SF_scalarField.h
/// @brief 与主流场同网格的通用标量场存储。

#include "SF_field.h"

#include <string>
#include <vector>

namespace SF {

/// @brief 绑定到 `Field` 网格尺寸的独立标量场。
///
/// 该类只负责存储和索引，不包含具体方程语义。后续 alpha、species、
/// temperature correction、turbulence scalar 和 passive tracer 都可以复用。
class ScalarField {
public:
    /// @brief 默认构造，不分配内存。
    ScalarField() = default;

    /// @brief 按主流场尺寸分配标量数组。
    /// @param field 参考主流场。
    /// @param name 标量名，用于诊断和 VTK 输出。
    /// @param fillValue 初始填充值。
    void setupLike(const Field& field,
                   const std::string& name,
                   double fillValue = 0.0);

    /// @brief 检查尺寸是否与给定主流场一致。
    /// @param field 参考主流场。
    /// @return 尺寸和存储均一致时返回 true。
    bool isCompatibleWith(const Field& field) const;

    /// @brief 查询是否尚未分配。
    bool empty() const { return values_.empty(); }

    /// @brief 标量名。
    const std::string& name() const { return name_; }

    /// @brief 重命名标量。
    /// @param name 新标量名。
    void setName(const std::string& name) { name_ = name; }

    /// @brief 含虚胞总维度。
    int MX() const { return mx_; }
    int MY() const { return my_; }
    int MZ() const { return mz_; }
    int NG() const { return ng_; }
    int TotalSize() const { return totalSize_; }

    /// @brief 计算平铺索引。
    /// @param i i方向索引。
    /// @param j j方向索引。
    /// @param k k方向索引。
    /// @return 0-based 平铺索引。
    int getIdx(int i, int j, int k) const {
        return (k * my_ + j) * mx_ + i;
    }

    /// @brief 可写访问。
    double& operator()(int i, int j, int k) {
        return values_[(size_t)getIdx(i, j, k)];
    }

    /// @brief 只读访问。
    double operator()(int i, int j, int k) const {
        return values_[(size_t)getIdx(i, j, k)];
    }

    /// @brief 直接访问底层数组。
    std::vector<double>& values() { return values_; }
    const std::vector<double>& values() const { return values_; }

private:
    std::string name_;
    int mx_ = 0;
    int my_ = 0;
    int mz_ = 0;
    int ng_ = 0;
    int totalSize_ = 0;
    std::vector<double> values_;
};

} // namespace SF
