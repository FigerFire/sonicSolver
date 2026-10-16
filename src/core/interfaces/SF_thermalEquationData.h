#pragma once
#include "core/field/SF_scalarField.h"
#include "core/state/SF_valueTypes.h"
namespace SF::FDM {
/// @brief 非拥有的温度物理存储与边界；数值 workspace 由 solver 创建。
struct ThermalEquationData {
    ScalarField* temperature=nullptr;
    std::vector<BCSetting<double>> boundaries;
};
}
