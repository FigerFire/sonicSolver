#pragma once

/// @file SF_reconstruction.h
/// @brief 物理边界定律与具体重构方法之间的统一调度接口。

#include "SF_field.h"
#include "law/SF_law.h"

namespace SF::Boundary::Reconstruction {

/// @brief 可选择的边界重构方法。
enum class Kind {
    Linear,
    Polynomial,
    MLS,
    ILW
};

/// @brief 一次边界重构选择。
struct Selection {
    Kind kind = Kind::Linear;
    int order = 0;
};

/// @brief 把现有显式 `[numerics].ILW` 配置转换为重构选择。
Selection fromILWSetting(bool enabled, int order);

/// @brief 对标量变量应用边界定律。
void applyScalar(Field& field, int i, int j, int k,
                 int axis, int variable, double prescribedValue,
                 Law::Kind law, Selection selection);

/// @brief 对 Vector3 变量应用边界定律。
void applyVector(Field& field, int i, int j, int k,
                 int axis, int variable, const Vector3& prescribedValue,
                 Law::Kind law, Selection selection);

} // namespace SF::Boundary::Reconstruction
