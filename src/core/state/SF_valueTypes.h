#pragma once

/// @file SF_valueTypes.h
/// @brief 跨模块传递的基础值类型，不依赖 parser 或求解器全局状态。

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>

#include "methods/math/frame/SF_vector.h"
#include "methods/math/SF_tensor.h"

namespace SF {

/// @brief 过渡兼容别名：纯数学权威已迁移到 SF::Math。
/// 新代码请直接使用 SF::Math::Vector3 / SF::Math::Tensor3 等。
using Vector3 = Math::Vector3;
using Tensor3 = Math::Tensor3;
using SymmTensor3 = Math::SymmTensor3;
using Math::cross;
using Math::dot;
using Math::norm;
using Math::normalize;

/// @brief 代数边界条件类型。
enum BCType { FIXED_VALUE, ZERO_GRADIENT, SYMMETRY, EMPTY };

/// @brief 求解域单元分类。
enum CellType { FLUID_CELL = 0, IBM_GHOST_CELL = 1, SOLID_CELL = 2 };

/// @brief 热边界条件类型。
enum class ThermalBCType {
    FixedTemperature,
    ZeroGradient,
    Adiabatic,
    HeatFlux,
    Empty
};

/// @brief 一个集合上的标量或向量条件。
template <typename T>
struct BCSetting {
    std::string name;
    BCType type;
    T value;
};

/// @brief 温度或热流边界条件。
struct ThermalBCSetting {
    std::string name;
    ThermalBCType type = ThermalBCType::ZeroGradient;
    double value = 0.0;
};

/// @brief 指定 zone 上的向量配置。
struct ZoneVectorSetting {
    std::string zone = "all";
    Vector3 value;
};

/// @brief 旋转坐标系配置。
struct RotatingSetting {
    std::string zone = "all";
    Vector3 center;
    Vector3 axis;
    double omega = 0.0;
    Vector3 velocity;
    bool hasVelocity = false;
};

/// @brief 壁面热源配置；heatFlux 为注入流体域的面热流。
struct WallHeatSetting {
    std::string patch = "all";
    std::string mode = "fixedHeatFlux";
    double heatFlux = 0.0;
    bool coupleEnergy = true;
    int rangeCoordinate = -1;
    double rangeMinimum = std::numeric_limits<double>::quiet_NaN();
    double rangeMaximum = std::numeric_limits<double>::quiet_NaN();
    double wallTemperature = std::numeric_limits<double>::quiet_NaN();
};

/// @brief 单流体五个守恒分量的兼容索引。
enum VarIdx { RHO = 0, RU = 1, RV = 2, RW = 3, E = 4, NVAR = 5 };

} // namespace SF
