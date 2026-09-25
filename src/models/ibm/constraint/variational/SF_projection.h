#pragma once

/// @file SF_projection.h
/// @brief 离散最小质量范数 IBM 投影及其驻值诊断。

#include "SF_valueTypes.h"
#include "SF_immersedConstraint.h"

#include <functional>
#include <vector>

namespace SF::IBM::Variational {

/// @brief 单个 Eulerian 体约束点的离散变分问题。
struct PointProjectionProblem {
    Vector3 predictedVelocity;
    Vector3 targetVelocity;
    double density = 0.0;
    double measure = 0.0;
    double dt = 0.0;
};

/// @brief `min 1/2 ||u-u*||_M^2, u=Us` 的驻值解。
struct PointProjectionSolution {
    Vector3 correctedVelocity;
    /// @brief 流体所受体力密度 lambda，单位 N/m3。
    Vector3 multiplier;
    double kineticIncrementFunctional = 0.0;
    double constraintResidual = 0.0;
    double stationarityResidual = 0.0;
};

/// @brief 精确求解单点质量矩阵 KKT，不使用 penalty 或隐藏降阶。
PointProjectionSolution solvePointProjection(
    const PointProjectionProblem& problem);

/// @brief 计算守恒动量校正对应的离散机械能增量。
double mechanicalEnergyIncrement(
    const Vector3& momentumIncrement,
    const Vector3& oldVelocity,
    const Vector3& newVelocity);

/// @brief 表面约束 `dt J M^-1 J^T M_L Lambda = U_s-Ju*` 的局部 owner-edge 视图。
struct SurfaceSchurProjectionProblem {
    const FDM::ImmersedSurfaceSystem* surface = nullptr;
    /// @brief 以 Eulerian local-cell 编号索引的 `1/(rho V)`。
    const std::vector<double>* inverseEulerianMass = nullptr;
    std::vector<Vector3> constraintError;
    double dt = 0.0;
    int maxIterations = 0;
    double relativeTolerance = 0.0;
    /// @brief 已复制的全局 marker 向量的唯一 owner 掩码。
    ///
    /// MPI 中 Schur image 与 rhs 在每 rank 都是完整副本；CG 内积只由一个
    /// marker owner 累加，随后 Runtime 执行全局 SUM，避免重复计数。
    /// 空数组保持串行与旧调用方的全 marker 语义。
    std::vector<unsigned char> ownedMarkers;
    /// @brief 可选的 Runtime SUM；未给出时即为严格串行算子。
    /// @details 数学核只调用此显式归约，不识别 MPI/rank/communicator。
    std::function<void(std::vector<double>&)> globalSum;
};

/// @brief 三个速度分量共享同一 Schur 算子的 matrix-free CG 解。
struct SurfaceSchurProjectionSolution {
    std::vector<Vector3> multiplier;
    int iterations = 0;
    double maximumRelativeResidual = 0.0;
};

/// @brief 精确考虑 marker 行重叠的表面 Schur 投影。
///
/// 用 `sqrt(M_L)` 对系统对称化；奇异、非正定或未收敛均 fail-fast，绝不
/// 静默改用 marker-local 对角近似。分布式时每个 rank 仅提供 owner edge 的
/// local contribution，`globalSum` 在每次 J/S 作用和内积处形成全局算子。
SurfaceSchurProjectionSolution solveSurfaceSchurProjection(
    const SurfaceSchurProjectionProblem& problem);

} // namespace SF::IBM::Variational
