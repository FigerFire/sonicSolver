#pragma once

/// @file SF_immersedStrategy.h
/// @brief IBM enforcement 的算法侧能力描述，不包含具体 Field/MPI 实现。

#include "SF_config.h"
#include "SF_immersedSystem.h"

#include <string>

namespace SF::ImmersedAlgorithm {

/// @brief 一个 enforcement strategy 对统一时间驱动器的静态契约。
struct StrategyContract {
    FDM::IBMEnforcement kind = FDM::IBMEnforcement::FractionalDLM;
    std::string name;
    bool requiresPredictedState = false;
    bool requiresPressureCorrection = false;
    bool producesMultiplier = false;
    bool implemented = false;
};

/// @brief 返回某一策略的数据契约；未知枚举直接失败。
StrategyContract contract(FDM::IBMEnforcement enforcement);

/// @brief 校验策略与顶层 flow algorithm 的数学阶段是否一致。
///
/// 这里只做阶段/能力检查，不负责调用 MPI，也不把策略替换成另一种策略。
void validateForAlgorithm(const FDM::IBMForcingConfig& forcing,
                          FDM::SolverAlgorithm algorithm);

/// @brief 以运行时统一方法选择校验 enforcement 所属的数学阶段。
///
/// Algorithm 只读取该中立 selection，不读取 IBMProperties，也不询问 legacy
/// constraint adapter 来猜测压力基/密度基路由。
void validateForAlgorithm(const FDM::ImmersedMethodSelection& selection,
                          FDM::SolverAlgorithm algorithm);

/// @brief 同时校验 selection 与模块实际声明的能力。
void validateForAlgorithm(const FDM::IImmersedSystem& system,
                          FDM::SolverAlgorithm algorithm);

} // namespace SF::ImmersedAlgorithm
