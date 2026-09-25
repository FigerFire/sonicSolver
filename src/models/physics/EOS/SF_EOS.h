/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.07.10-----------*/

#pragma once

/// @file SF_EOS.h
/// @brief 相/组分状态方程统一接口。

#include <memory>
#include <stdexcept>
#include <string>

namespace SF::Physics::EOS {

/// @brief 单个相—组分的热力学状态。
struct State {
    double density = 0.0;
    double internalEnergy = 0.0;
    double pressure = 0.0;
    double temperature = 0.0;
    double soundSpeed = 0.0;
};

/// @brief EOS 只负责热力学映射，不读写 Field。
class Model {
public:
    virtual ~Model() = default;
    virtual std::string type() const = 0;
    virtual State fromDensityEnergy(double density,
                                    double internalEnergy) const = 0;
    /// @brief 固定密度和压力反演比内能，用于压力边界与冷启动。
    virtual State fromDensityPressure(double, double) const {
        throw std::runtime_error(
            "EOS does not implement fixed-density pressure inversion.");
    }
    /// @brief 固定密度和温度反演比内能，用于温度边界与冷启动。
    virtual State fromDensityTemperature(double, double) const {
        throw std::runtime_error(
            "EOS does not implement fixed-density temperature inversion.");
    }
    virtual State fromPressureTemperature(double pressure,
                                          double temperature) const = 0;
    virtual double referenceEnergy() const = 0;
};

using ModelPtr = std::shared_ptr<const Model>;

} // namespace SF::Physics::EOS
