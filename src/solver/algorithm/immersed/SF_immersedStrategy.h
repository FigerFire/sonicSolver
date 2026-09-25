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

/// @brief 校验顺序约束投影 operation 的方法与能力。
void validateProjectionProvider(const FDM::IImmersedSystem& system);

/// @brief 校验 monolithic KKT operation 的方法与能力。
void validateMonolithicProvider(const FDM::IImmersedSystem& system);

} // namespace SF::ImmersedAlgorithm
