/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.06.07-----------*/

#pragma once

/// @file SF_thermodynamicClosure.h
/// @brief 标量热力学边界适配：压力输入转化为守恒总能量闭合。

#include "SF_field.h"
#include "SF_idealGasDefaults.h"

#include <vector>

namespace SF {
namespace Boundary {

/// @brief 从守恒量计算压力。
/// @param field 守恒量场。
/// @param i i索引。
/// @param j j索引。
/// @param k k索引。
/// @return 由当前守恒状态得到的压力。
double pressureAt(const Field& field, int i, int j, int k);

/// @brief 从守恒量计算温度。
/// @param field 守恒量场。
/// @param i i索引。
/// @param j j索引。
/// @param k k索引。
/// @return 理想气体温度。
double temperatureAt(const Field& field, int i, int j, int k);

/// @brief 用指定压力和当前密度/动量重建守恒总能量。
/// @param field 守恒量场。
/// @param i i索引。
/// @param j j索引。
/// @param k k索引。
/// @param pressure 边界压力。
void setEnergyFromPressure(Field& field, int i, int j, int k,
                           double pressure);

/// @brief 用指定温度和当前密度/动量重建守恒总能量。
/// @param field 守恒量场。
/// @param i i索引。
/// @param j j索引。
/// @param k k索引。
/// @param temperature 边界温度。
void setEnergyFromTemperature(Field& field, int i, int j, int k,
                              double temperature);

/// @brief 按压力边界条件闭合，并把结果回写到守恒总能量 E。
/// @param field 守恒量场。
/// @param pSettings IO读取的压力边界条件。
/// @param useILW true时使用普通物理边界ILW闭合。
void updateEnergyFromPressure(Field& field,
                              const std::vector<BCSetting<double>>& pSettings,
                              bool useILW = false,
                              int ilwOrder = 0);

/// @brief 按温度/热流边界条件闭合，并把结果回写到守恒总能量 E。
/// @param field 守恒量场。
/// @param thermalSettings IO读取的温度/热流边界条件。
/// @param dynamicViscosity 常动力粘度。
/// @param prandtl Prandtl 数。
/// @param useILW true时当前会 fail-fast，避免隐藏降阶热边界。
void updateEnergyFromThermalBoundary(
    Field& field,
    const std::vector<ThermalBCSetting>& thermalSettings,
    double dynamicViscosity,
    double prandtl,
    bool useILW = false);

} // namespace Boundary
} // namespace SF
