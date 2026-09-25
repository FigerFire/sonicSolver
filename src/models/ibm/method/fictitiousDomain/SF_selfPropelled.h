#pragma once

/// @file SF_selfPropelled.h
/// @brief DFM 自推进算法共享的广义刚体速度离散求解器。

#include "SF_configTypes.h"

#include <array>
#include <vector>

namespace SF::IBM::Forcing::SelfPropelled {

/// @brief 一个虚拟域点对平动/转动零合力约束的贡献。
struct Sample {
    Vector3 predictedVelocity;
    Vector3 deformationVelocity;
    Vector3 relativePosition;
    double phaseMass = 0.0;
};

/// @brief 不含外载荷的局部广义质量系统贡献；可由 Runtime 逐项 SUM。
struct GeneralizedSystem {
    std::array<double,9> matrix{};
    Vector3 rhs;
};

/// @brief 从本 rank 唯一 Eulerian body 点装配局部广义质量贡献。
GeneralizedSystem assembleGeneralizedSystem(
    const std::vector<Sample>& samples,
    FDM::IBMRigidMotionMode mode);

/// @brief 从已经全局归并的广义系统恢复刚体速度。
Vector3 solveGeneralizedSystem(
    const GeneralizedSystem& system,
    FDM::IBMRigidMotionMode mode,
    const Vector3& externalForce,
    const Vector3& externalTorque,
    double dt);

/// @brief 由 Bhalla DFM 的零总乘子条件解三个广义速度分量。
/// @return `motivation` 时返回平动速度，`rotate` 时返回角速度。
Vector3 solveGeneralizedVelocity(
    const std::vector<Sample>& samples,
    FDM::IBMRigidMotionMode mode,
    const Vector3& externalForce,
    const Vector3& externalTorque,
    double dt);

/// @brief 把广义速度与局部形变速度组合成实体点目标速度。
Vector3 targetVelocity(
    FDM::IBMRigidMotionMode mode,
    const Vector3& generalizedVelocity,
    const Vector3& relativePosition,
    const Vector3& deformationVelocity);

} // namespace SF::IBM::Forcing::SelfPropelled
