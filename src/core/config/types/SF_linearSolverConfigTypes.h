#pragma once
/// @file SF_linearSolverConfigTypes.h
/// @brief 大型稀疏线性求解后端的强类型配置值对象。

#include <cstdint>
#include <cmath>
#include <stdexcept>
#include <string>

namespace SF::FDM {
/// @brief 大型稀疏系统后端。
enum class LinearBackend { Hypre };
/// @brief HYPRE Krylov 方法。
enum class KrylovMethod { FlexGMRES, PCG };
/// @brief HYPRE 预条件器族。
enum class LinearPreconditioner { BoomerAMG, ILU, KKTBlockSchur };
/// @brief 送入线性后端前的显式代数平衡方式。
enum class LinearEquilibration { None, RowMax };

/// @brief 单个方程族的大型线性求解配置。
struct LinearSolverConfig {
    LinearBackend backend = LinearBackend::Hypre;
    KrylovMethod method = KrylovMethod::FlexGMRES;
    LinearPreconditioner preconditioner = LinearPreconditioner::BoomerAMG;
    LinearEquilibration equilibration = LinearEquilibration::None;
    int maxIterations = 500;
    int krylovDimension = 50;
    double relativeTolerance = 1.0e-8;
    double absoluteTolerance = 1.0e-12;
    int structureRebuildInterval = 0;
    int preconditionerRefreshInterval = 10;
    int iluType = 0;
    int iluLevelOfFill = 0;
    /// @brief BoomerAMG 显式参数；也用于 KKT Schur 压力子块。
    int amgCoarsenType = 10;
    int amgRelaxType = 6;
    int amgSweeps = 1;
    int amgMaxLevels = 25;
    double amgStrongThreshold = 0.25;
    /// @brief 串行稠密约束 Schur 的容量限制，超限必须报错。
    int schurDenseLimit = 256;
};

/// @brief 校验线性求解配置，不允许无效容差或隐式算法替换。
///
/// 这是 typed 值对象的自校验，只依赖本文件声明的枚举与字段，因此属于 core，
/// 不属于 application 解析层。
/// @param config 待校验的线性求解配置。
/// @param context 诊断上下文（方程族名）。
inline void validateLinearSolverConfig(const LinearSolverConfig& config,
                                       const std::string& context) {
    if (config.maxIterations <= 0 || config.krylovDimension <= 0
        || !std::isfinite(config.relativeTolerance)
        || config.relativeTolerance <= 0.0
        || !std::isfinite(config.absoluteTolerance)
        || config.absoluteTolerance < 0.0
        || config.structureRebuildInterval < 0
        || config.preconditionerRefreshInterval < 0
        || config.iluLevelOfFill < 0
        || config.amgSweeps <= 0 || config.amgMaxLevels <= 0
        || config.amgCoarsenType < 0 || config.amgRelaxType < 0
        || !std::isfinite(config.amgStrongThreshold)
        || config.amgStrongThreshold <= 0 || config.amgStrongThreshold > 1
        || config.schurDenseLimit <= 0) {
        throw std::invalid_argument(
            context + ": invalid Krylov iteration/tolerance controls.");
    }
    if (config.method == KrylovMethod::PCG
        && config.preconditioner != LinearPreconditioner::BoomerAMG) {
        throw std::invalid_argument(
            context
            + ": PCG requires the symmetric boomerAMG preconditioner.");
    }
    if (config.method == KrylovMethod::PCG
        && config.equilibration != LinearEquilibration::None) {
        throw std::invalid_argument(
            context
            + ": rowMax scaling is not symmetric; use flexGMRES.");
    }
    if (config.preconditioner == LinearPreconditioner::ILU
        && config.iluType != 0 && config.iluType != 1
        && config.iluType != 10 && config.iluType != 11
        && config.iluType != 20 && config.iluType != 21
        && config.iluType != 30 && config.iluType != 31
        && config.iluType != 40 && config.iluType != 41
        && config.iluType != 50) {
        throw std::invalid_argument(
            context + ": unsupported HYPRE ILU type.");
    }
}
} // namespace SF::FDM
