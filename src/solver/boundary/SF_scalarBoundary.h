#pragma once
/// @file SF_scalarBoundary.h
/// @brief 节点标量边界：fixedValue 镜像，多层 zeroGradient 常值延拓，EMPTY 两平面复制。
#include "core/field/SF_field.h"
#include "core/field/SF_scalarField.h"
#include "SF_config.h"
namespace SF::Boundary {
/// @brief 写入原标量边界和 ghost；不执行闭合、求解或 halo 通信。
void applyNodeScalarBoundary(const Field& geometry,
    const std::vector<BCSetting<double>>& settings,ScalarField& scalar);
/// @brief 既有single-fluid输运BC；保持多轴corner和low-side ghost写入顺序。
void applyTransportScalarBoundary(const Field&,const std::vector<BCSetting<double>>&,std::vector<double>&);
}
