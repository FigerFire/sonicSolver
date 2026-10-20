/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.07.07-----------*/

#pragma once

/// @file SF_phaseChange.h
/// @brief 相变模型调度入口。

#include "SF_field.h"
#include "methods/numerics/structured/SF_structured.h"
#include "SF_scalarField.h"
#include "SF_phaseProperties.h"

#include <limits>
#include <string>
#include <vector>

namespace SF {
namespace Physics {
namespace PhaseChange {

/// @brief 单速度 mixture 相变模型的局部上下文。
struct ModelContext {
    const Field& field;
    const ScalarField& alpha;
    const ScalarField& temperature;
    const Multiphase::MultiPhaseConfig& config;
    const Multiphase::PhaseProperties& alphaPhase;
    const Multiphase::PhaseProperties& otherPhase;
    const Multiphase::PhaseProperties& liquidPhase;
    bool alphaIsLiquid = true;
    double dt = 0.0;
    double localPressure = std::numeric_limits<double>::quiet_NaN();
    double localLiquidSpeed = std::numeric_limits<double>::quiet_NaN();
    double localWallDistance = std::numeric_limits<double>::quiet_NaN();
    double localFrictionVelocity = std::numeric_limits<double>::quiet_NaN();
};

/// @brief 相变源项在每个时间步内的聚合诊断量。
struct PhaseChangeDiagnostics {
    int totalFluidCells = 0;
    int activeCells = 0;           ///< mdot != 0 的格子数
    int limitedCells = 0;          ///< 被 alpha/热量保守限制的格子数
    int failedMassCells = 0;       ///< 限制后仍无可用相质量的格子数
    int failedHeatCells = 0;       ///< 限制后仍无可用热量的格子数
    int failedPostStateCells = 0;  ///< 源项施加后物理状态无效的格子数
    double maxRateReduction = 0.0; ///< 最大限制比例 |rate-actual|/|rate|
    double minAlpha = 2.0;
    double maxAlpha = -1.0;
    double minTemperature = 1.0e300;
    double maxTemperature = -1.0e300;
    double minMdot = 1.0e300;
    double maxMdot = -1.0e300;
    double minEnergySource = 1.0e300;
    double maxEnergySource = -1.0e300;
    int firstBadCellIndex = -1;
};

/// @brief 相变定律给出的候选交换率；不负责守恒变量映射或相质量限制。
struct PhaseChangeRateResult {
    std::vector<double> mdot;           ///< 候选质量传递速率 (kg/(m^3*s))。
    PhaseChangeDiagnostics diagnostics; ///< 聚合诊断信息。
};

/// @brief 模型名称规整。
std::string normalizeModel(std::string model);

/// @brief 读取模型专用系数。
double coefficient(const Multiphase::PhaseChangeOptions& config,
                   const std::string& key,
                   double fallback = 0.0);

/// @brief 读取模型专用字符串选项。
std::string selection(const Multiphase::PhaseChangeOptions& config,
                      const std::string& key,
                      const std::string& fallback = "");

/// @brief 校验 phaseChange 配置。
void validateConfig(const Multiphase::MultiPhaseConfig& config,
                    const std::string& context);

/// @brief 计算局部相变候选速率。
///
/// 质量源符号统一为 `mdot > 0` 表示 liquid -> vapour。该入口不按 alpha、
/// 配置密度或 dt 限制速率，也不生成 alpha/rhoE 源；FluidStateModel 的 transfer
/// assembler 必须使用当前热力学状态完成质量耗尽限制和守恒映射。
///
/// @param field             守恒量场（只读）。
/// @param alpha             alpha 标量场（只读）。
/// @param temperature       温度标量场（只读）。
/// @param config            多相配置。
/// @param dt                时间步长。
/// @return 候选质量交换率和诊断量。
PhaseChangeRateResult computeRates(const Field& field,
                                   const ScalarField& alpha,
                                   const ScalarField& temperature,
                                   const Multiphase::MultiPhaseConfig& config,
                                   double dt);

/// @brief 仅在时间层提交后输出已缓存的相变诊断。
void reportDiagnostics(const PhaseChangeDiagnostics& diagnostics);

// ---- 各模型的质量速率计算入口 ----
void addLeeRates(const ModelContext& ctx, std::vector<double>& mdot);
void addStefanRates(const ModelContext& ctx, std::vector<double>& mdot);
void addHKSRates(const ModelContext& ctx, std::vector<double>& mdot);
void addEnthalpyPorosityRates(const ModelContext& ctx,
                              std::vector<double>& mdot);
void addSchnerrSauerRates(const ModelContext& ctx, std::vector<double>& mdot);
void addZGBRates(const ModelContext& ctx, std::vector<double>& mdot);
void addPhaseFieldRates(const ModelContext& ctx, std::vector<double>& mdot);
void addSaturationPropertyRates(const ModelContext& ctx,
                                std::vector<double>& mdot);

} // namespace PhaseChange
} // namespace Physics
} // namespace SF
