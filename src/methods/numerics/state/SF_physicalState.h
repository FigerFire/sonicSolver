/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.07.08-----------*/

#pragma once

/// @file SF_physicalState.h
/// @brief 守恒状态热力学闭合诊断。

#include "SF_field.h"
#include "methods/numerics/structured/SF_structured.h"
#include "SF_utility.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

namespace SF {
namespace Numerics {

/// @brief 一个物理量最小值及其网格位置。
struct StateMinimum {
    double value = std::numeric_limits<double>::max(); ///< 最小值。
    int i = 0; ///< i方向索引。
    int j = 0; ///< j方向索引。
    int k = 0; ///< k方向索引。
};

/// @brief 一次守恒状态闭合检查的聚合结果。
struct PhysicalStateSummary {
    int checkedCells = 0;       ///< 被检查的真实流体点数量。
    StateMinimum rho;           ///< 最小密度。
    StateMinimum pressure;      ///< 最小压力。
};

/// @brief 记录一个最小值。
/// @param minimum 当前最小值记录。
/// @param value 候选值。
/// @param i i方向索引。
/// @param j j方向索引。
/// @param k k方向索引。
inline void updateMinimum(StateMinimum& minimum,
                          double value,
                          int i,
                          int j,
                          int k) {
    if (value < minimum.value) {
        minimum.value = value;
        minimum.i = i;
        minimum.j = j;
        minimum.k = k;
    }
}

/// @brief 输出非法守恒状态并终止。
/// @param context 调用位置或RK stage名称。
/// @param field 守恒量场。
/// @param i i方向索引。
/// @param j j方向索引。
/// @param k k方向索引。
/// @param pressureValue 已计算压力；无法计算时为NaN。
inline void fatalInvalidPhysicalState(const char* context,
                                      const Field& field,
                                      int i,
                                      int j,
                                      int k,
                                      double pressureValue) {
    std::cerr << "[SF FATAL] " << context
              << ": invalid conservative state after closure check at ("
              << i << "," << j << "," << k << ")"
              << " rho=" << field(i, j, k, RHO)
              << " ru=" << field(i, j, k, RU)
              << " rv=" << field(i, j, k, RV)
              << " rw=" << field(i, j, k, RW)
              << " E=" << field(i, j, k, E)
              << " p=" << pressureValue
              << std::endl;
    std::exit(1);
}

/// @brief 检查守恒状态并返回聚合诊断。
///
/// 本函数只检查密度基守恒量和压力，不在 Euler 数值层猜测温度。
/// 温度必须由所属热力学模型闭合和诊断。
///
/// @param field 守恒量场。
/// @param context 调用位置或RK stage名称。
/// @return 聚合后的最小值统计。
inline PhysicalStateSummary validatePhysicalState(const Field& field,
                                                  const char* context) {
    PhysicalStateSummary summary;
    Math::forFluidInterior(field, [&](int i, int j, int k) {
        const double rho = field(i, j, k, RHO);
        const double ru = field(i, j, k, RU);
        const double rv = field(i, j, k, RV);
        const double rw = field(i, j, k, RW);
        const double e = field(i, j, k, E);
        const bool badConservative =
            !finite(rho) || !finite(ru) || !finite(rv)
            || !finite(rw) || !finite(e) || rho <= 0.0;
        const double p = badConservative
            ? std::numeric_limits<double>::quiet_NaN()
            : pressure(rho, ru, rv, rw, e);
        if (badConservative || !finite(p) || p <= 0.0) {
            fatalInvalidPhysicalState(context, field, i, j, k, p);
        }

        ++summary.checkedCells;
        updateMinimum(summary.rho, rho, i, j, k);
        updateMinimum(summary.pressure, p, i, j, k);
    });

    if (summary.checkedCells == 0) {
        std::cerr << "[SF FATAL] " << context
                  << ": physical state closure found no solved fluid cells."
                  << std::endl;
        std::exit(1);
    }
    return summary;
}

/// @brief 生成闭合诊断摘要文本。
/// @param context 调用位置或RK stage名称。
/// @param summary 聚合统计。
/// @return 可写入日志的单行摘要。
inline std::string formatPhysicalStateSummary(
        const char* context,
        const PhysicalStateSummary& summary) {
    std::ostringstream os;
    os << context
       << ": checked=" << summary.checkedCells
       << ", min(rho)=" << summary.rho.value << "@("
       << summary.rho.i << "," << summary.rho.j << ","
       << summary.rho.k << ")"
       << ", min(p)=" << summary.pressure.value << "@("
       << summary.pressure.i << "," << summary.pressure.j << ","
       << summary.pressure.k << ")";
    return os.str();
}

} // namespace Numerics
} // namespace SF
