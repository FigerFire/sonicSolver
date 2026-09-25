/// @file SF_state.h
/// @brief Level Set 输运、重初始化或界面几何模型实现。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#pragma once

#include "SF_field.h"
#include "SF_phaseProperties.h"
#include "SF_scalarField.h"

#include <functional>
#include <string>
#include <vector>

namespace SF {
namespace Physics {
namespace Multiphase {

/// @brief Level Set 标量场、界面几何和相态物性缓存。
///
/// phi > 0 表示液相, phi < 0 表示气相。该类只保存多相变量,
/// 不直接修改守恒量 Field。
class LevelSetField {
public:
    /// @brief 按 Field 尺寸分配 Level Set 存储。
    /// @param field 参考流场。
    /// @param fillValue 初始填充值。
    void setupLike(const Field& field, double fillValue = 0.0);

    /// @brief 用 multiPhase 配置中的默认相和 set 相赋值初始化符号场。
    /// @param field 已加载网格和 set 的流场。
    /// @param config 多相配置。
    void initializeFromSets(const Field& field,
                            const MultiPhaseConfig& config);

    /// @brief 按 0/phi 中的 boundaryField 条件更新 Level Set 边界值。
    /// @param field 已加载网格和 set 的流场。
    /// @param config 多相配置。
    void applyBoundaryConditions(const Field& field,
                                 const MultiPhaseConfig& config);

    /// @brief 计算界面法向和曲率缓存。
    /// @param field 参考网格。
    /// @param gradientTolerance 用户输入的最小有效梯度。
    void computeGeometry(const Field& field, double gradientTolerance);

    /// @brief 解重初始化方程, 将 phi 推向 signed-distance 场。
    /// @param field 参考网格。
    /// @param pseudoSteps 伪时间步数, 必须大于0。
    /// @param pseudoTimeStep 伪时间步长, 必须为正。
    /// @param order WENO阶数, 仅支持3/5/7。
    /// @param wenoEpsilon,wenoPower HJ-WENO 非线性权重参数。
    /// @param signSmoothingFactor 符号函数平滑宽度相对网格尺度倍率。
    void reinitialize(const Field& field,
                      int pseudoSteps,
                      double pseudoTimeStep,
                      int order,
                      double wenoEpsilon,
                      double wenoPower,
                      double signSmoothingFactor,
                      const std::function<void()>& prepareStage = {});

    /// @brief 根据phi更新相态密度和黏度缓存。
    /// @param field 参考网格。
    /// @param config 多相配置。
    void updateMaterialProperties(const Field& field,
                                  const MultiPhaseConfig& config);

    /// @brief 检查是否已经按给定 Field 完成分配。
    bool isCompatibleWith(const Field& field) const;

    /// @brief 含虚胞的总 X/Y/Z 网格数。
    int MX() const { return mx_; }
    int MY() const { return my_; }
    int MZ() const { return mz_; }
    int NG() const { return ng_; }
    int TotalSize() const { return totalSize_; }

    /// @brief 计算 Level Set 平铺索引。
    int getIdx(int i, int j, int k) const { return (k * my_ + j) * mx_ + i; }

    /// @brief 访问 phi。
    double& phi(int i, int j, int k) { return phi_(i, j, k); }
    double phi(int i, int j, int k) const { return phi_(i, j, k); }

    /// @brief 访问缓存法向。
    const Vector3& normal(int i, int j, int k) const {
        return normals_[getIdx(i, j, k)];
    }

    /// @brief 访问缓存曲率。
    double curvature(int i, int j, int k) const {
        return curvature_[getIdx(i, j, k)];
    }

    /// @brief 查询该点是否处于符号变化界面邻域。
    bool isInterfaceCell(int i, int j, int k) const {
        return interfaceMask_[getIdx(i, j, k)] != 0;
    }

    /// @brief 查询该点是否在窄带内。
    bool inNarrowBand(int i, int j, int k) const {
        return narrowBandMask_[getIdx(i, j, k)] != 0;
    }

    /// @brief 查询相态物性缓存是否已更新。
    bool hasMaterialProperties() const { return hasMaterialProperties_; }

    /// @brief 访问混合密度缓存。
    double density(int i, int j, int k) const {
        return density_[getIdx(i, j, k)];
    }

    /// @brief 访问混合动力黏度缓存。
    double viscosity(int i, int j, int k) const {
        return viscosity_[getIdx(i, j, k)];
    }

    /// @brief 直接访问 phi 存储, 供数值核批量更新。
    std::vector<double>& values() { return phi_.values(); }
    const std::vector<double>& values() const { return phi_.values(); }

    /// @brief 返回统一状态注册表可消费的 Level Set 主标量。
    ScalarField& scalar() { return phi_; }
    const ScalarField& scalar() const { return phi_; }

    /// @brief 直接访问法向缓存, 供界面几何核写入。
    std::vector<Vector3>& normals() { return normals_; }
    const std::vector<Vector3>& normals() const { return normals_; }

    /// @brief 直接访问曲率缓存, 供界面几何核写入。
    std::vector<double>& curvatures() { return curvature_; }
    const std::vector<double>& curvatures() const { return curvature_; }

    /// @brief 直接访问界面掩码, 供界面几何核写入。
    std::vector<unsigned char>& interfaceMask() { return interfaceMask_; }
    const std::vector<unsigned char>& interfaceMask() const { return interfaceMask_; }

    /// @brief 直接访问窄带掩码, 供窄带模块写入。
    std::vector<unsigned char>& narrowBandMask() { return narrowBandMask_; }
    const std::vector<unsigned char>& narrowBandMask() const { return narrowBandMask_; }

    /// @brief 直接访问混合密度缓存。
    std::vector<double>& densities() { return density_; }
    const std::vector<double>& densities() const { return density_; }

    /// @brief 直接访问混合黏度缓存。
    std::vector<double>& viscosities() { return viscosity_; }
    const std::vector<double>& viscosities() const { return viscosity_; }

    /// @brief 标记物性缓存状态。
    /// @param enabled true表示density/viscosity缓存已与当前phi一致。
    void setMaterialPropertiesReady(bool enabled) {
        hasMaterialProperties_ = enabled;
    }

private:
    int mx_ = 0;
    int my_ = 0;
    int mz_ = 0;
    int ng_ = 0;
    int totalSize_ = 0;
    ScalarField phi_;
    std::vector<Vector3> normals_;
    std::vector<double> curvature_;
    std::vector<unsigned char> interfaceMask_;
    std::vector<unsigned char> narrowBandMask_;
    std::vector<double> density_;
    std::vector<double> viscosity_;
    bool hasMaterialProperties_ = false;
};

} // namespace Multiphase
} // namespace Physics
} // namespace SF
