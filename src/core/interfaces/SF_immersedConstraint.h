#pragma once

/// @file SF_immersedConstraint.h
/// @brief 浸没体约束校正与拉格朗日乘子诊断的求解器中立接口。

#include "core/field/SF_field.h"
#include "core/model/SF_globalEntityId.h"
#include "SF_valueTypes.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace SF::FDM {

class IImmersedSystem;
class IExecutionRuntime;

/// @brief 一次浸没约束校正的守恒诊断。
struct ImmersedConstraintResult {
    bool performed = false;
    std::size_t constrainedCells = 0;
    double maximumVelocityResidual = 0.0;
    /// @brief 离散驻值方程 `M(u-u*)-C^T lambda` 的最大残差。
    double maximumStationarityResidual = 0.0;
    /// @brief 本次约束校正的离散质量范数泛函值。
    double kineticIncrementFunctional = 0.0;
    Vector3 forceOnBody;
    Vector3 torqueOnBody;
    double fluidMechanicalPower = 0.0;
    std::string detail;
};

/// @brief 一个 Eulerian 体约束自由度的 canonical 单 patch 描述。
///
/// `localCell` 只在当前尚未分布式化的 Field 内有效；它不是 GlobalDofId，
/// 更不是 HYPRE row。分布式实现必须由 Runtime/mesh ownership 提供全局身份，
/// 不能由 patch 局部索引猜测。
struct ImmersedBodyPoint {
    int localCell = -1;
    Vector3 position;
    Vector3 relativePosition;
    double dualVolume = 0.0;
};

/// @brief body mask 约束的几何/拓扑数据，不包含流体状态或乘子。
struct ImmersedBodySystem {
    std::vector<ImmersedBodyPoint> points;
    Vector3 centerOfMass;
};

/// @brief 一个表面点对 Eulerian 速度自由度的插值权重。
struct ImmersedInterpolationWeight {
    int cell = -1;
    double value = 0.0;
    /// @brief canonical Eulerian entity identity；串行时由唯一网格点编号提供。
    std::int64_t globalEulerianDofId = -1;
    /// @brief 该 Eulerian 点的对偶体积，单位 m3。
    double dualVolume = 0.0;
};

/// @brief 一个 Lagrangian 表面积分点及其 J 行。
struct ImmersedSurfacePoint {
    /// @brief 稳定约束数学身份；不是 HYPRE row 或 marker 局部数组下标。
    GlobalConstraintDofId globalConstraintId;
    Vector3 position;
    Vector3 relativePosition;
    /// @brief prescribed solid 时的该积分点速度。
    Vector3 prescribedVelocity;
    /// @brief 把 `q=(U_s,omega_s)` 六个广义速度映射到点速度的列向量 G。
    std::array<Vector3,6> solidVelocityBasis;
    double measure = 0.0;
    std::vector<ImmersedInterpolationWeight> interpolation;
};

/// @brief KKT 可直接消费的固体方程 `A_s q = r_s`，数组槽位保持三平动三转动。
///
/// 该数据已经包含本时间步的质量/惯量、外载荷和陀螺项。线性系统装配器
/// 不需要知道 q 表示平动还是转动，也不读取 IBM 配置。
struct ImmersedSolidEquation {
    bool solveGeneralizedVelocity = false;
    int generalizedDofs = 0;
    /// @brief 固体模型允许的物理分量；空间维度约束再与此集合取交集。
    std::vector<int> activeComponents;
    std::array<double,36> lhs{};
    std::array<double,6> rhs{};
    std::array<double,6> initialGuess{};
};

/// @brief 压力算法装配单体 KKT 所需的只读表面和刚体状态。
struct ImmersedSurfaceSystem {
    std::vector<ImmersedSurfacePoint> points;
    ImmersedSolidEquation solidEquation;
    double constraintTolerance = 0.0;
    /// @brief 大于零时启用 `gamma/2 ||J du-r||^2_ML` 增广项。
    double augmentationCoefficient = 0.0;
};

/// @brief 单体 KKT 解回传给 IBM 模块的表面乘子与刚体速度。
struct ImmersedKKTState {
    std::vector<Vector3> surfaceMultiplier;
    bool hasGeneralizedVelocity = false;
    std::array<double,6> generalizedVelocity{};
    /// @brief KKT 动量驻值行 `M du + G^T p - J^T lambda = 0` 的最大残差。
    double maximumStationarityResidual = 0.0;
    /// @brief 速度增量对应的离散质量范数泛函值。
    double kineticIncrementFunctional = 0.0;
};

/// @brief predictor 与 solid kinematics 之间的代数约束接口。
///
/// 该接口不属于 ghost 边界管线，也不把乘子伪装成预先给定的体源。
/// 实现从离散驻值条件中求出满足 `J u = U_s` 所需的乘子，并一致更新状态。
class IImmersedConstraint {
public:
    virtual ~IImmersedConstraint() = default;

    /// @brief 返回创建该约束适配器的选择服务；未知时返回 nullptr。
    ///
    /// 当前 application adapter 返回同一个 `IImmersedSystem` 实例，算法可
    /// 据此拒绝把一个 manager 的几何约束和另一个 manager 的选择轴拼接。
    /// 老的外部实现可以保持 nullptr，不会因此获得新的数值 fallback。
    virtual const IImmersedSystem* systemProvider() const { return nullptr; }

    /// @brief 注入 Algorithm 当前使用的抽象 Runtime；实现不得保留 MPI/communicator。
    virtual void setExecutionRuntime(IExecutionRuntime*) {}

    /// @brief 在目标时间层对预测状态施加浸没约束。
    /// @param fields 本 rank 当前推进的 canonical 流体状态。
    /// @param targetTime 校正完成后的物理时间，单位 s。
    /// @param dt 本次时间增量，单位 s。
    virtual ImmersedConstraintResult projectPredictedState(
        const std::vector<Field*>& fields,
        double targetTime,
        double dt) = 0;

    /// @brief 构造当前时间层的表面 J/S 与刚体状态。
    virtual const ImmersedSurfaceSystem& prepareMonolithicSystem(
        Field&, double, double) {
        throw std::runtime_error(
            "Immersed constraint does not provide a monolithic KKT system.");
    }

    /// @brief 接收单体解、更新乘子诊断及刚体状态。
    virtual ImmersedConstraintResult acceptMonolithicSolution(
        Field&, double, double, const ImmersedKKTState&) {
        throw std::runtime_error(
            "Immersed constraint cannot accept a monolithic KKT solution.");
    }
};

} // namespace SF::FDM
