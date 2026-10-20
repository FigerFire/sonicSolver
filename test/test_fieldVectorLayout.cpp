/// @file test_fieldVectorLayout.cpp
/// @brief Field 变量组读写 contract：Vector3 的布局必须由传入的 v 决定。
///
/// 该测试保护一个真实的 correctness 缺陷：`getVal<Vector3>` / `setVal<Vector3>`
/// 曾经把分量固定写成 1/2/3（RU/RV/RW）。任何把 vector 组放在其他 offset 的
/// state model（multiphase、Eulerian、带 auxiliary 未知量的布局）都会静默读写
/// 错位。这里只检查布局 contract，不涉及任何数值公式。

#include "SF_field.h"

#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

/// @brief 对一个 storage layout 验证非 1/2/3 的 vector 组读写。
void checkVectorGroup(SF::FieldStorageLayout layout) {
    const int variableCount = 7;
    const int group = 2;  ///< 组内分量 2/3/4，刻意避开 legacy 1/2/3。
    SF::Field field;
    field.setStorageLayout(layout);
    field.setup(2, 2, 1, 1, variableCount);

    const int i = 1, j = 1, k = 1;
    for (int v = 0; v < variableCount; ++v) {
        field(i, j, k, v) = 100.0 + static_cast<double>(v);
    }

    const SF::Vector3 read = field.getVal<SF::Vector3>(i, j, k, group);
    require(read.x == field(i, j, k, group),
            "getVal<Vector3> must read component v as x.");
    require(read.y == field(i, j, k, group + 1),
            "getVal<Vector3> must read component v+1 as y.");
    require(read.z == field(i, j, k, group + 2),
            "getVal<Vector3> must read component v+2 as z.");

    field.setVal<SF::Vector3>(i, j, k, group, SF::Vector3(1.0, 2.0, 3.0));
    require(field(i, j, k, group) == 1.0, "setVal<Vector3> must write x at v.");
    require(field(i, j, k, group + 1) == 2.0, "setVal<Vector3> must write y at v+1.");
    require(field(i, j, k, group + 2) == 3.0, "setVal<Vector3> must write z at v+2.");
    for (int v = 0; v < variableCount; ++v) {
        if (v >= group && v <= group + 2) continue;
        require(field(i, j, k, v) == 100.0 + static_cast<double>(v),
                "setVal<Vector3> must not touch other variable groups.");
    }

    // 标量组仍按 v 读写。
    field.setVal<double>(i, j, k, 0, 7.5);
    require(field.getVal<double>(i, j, k, 0) == 7.5,
            "scalar getVal/setVal must use v directly.");
}

/// @brief 验证 legacy 布局（组起点 1）与新布局给出相同的组语义。
void checkLegacyMomentumOffset() {
    SF::Field field;
    field.setup(2, 2, 1, 1, 5);
    const int i = 1, j = 1, k = 1;
    field.setVal<SF::Vector3>(i, j, k, 1, SF::Vector3(0.5, -1.5, 2.5));
    require(field.getVal<SF::Vector3>(i, j, k, 1).x == 0.5
                && field.getVal<SF::Vector3>(i, j, k, 1).y == -1.5
                && field.getVal<SF::Vector3>(i, j, k, 1).z == 2.5,
            "legacy momentum group at v=1 must round-trip through getVal/setVal.");
    require(field(i, j, k, 0) == 0.0 && field(i, j, k, 4) == 0.0,
            "momentum group writes must not leak into density/energy.");
}

} // namespace

int main() {
    try {
        checkVectorGroup(SF::FieldStorageLayout::AoS);
        checkVectorGroup(SF::FieldStorageLayout::SoA);
        checkLegacyMomentumOffset();
    } catch (const std::exception& error) {
        std::cerr << "[test_fieldVectorLayout] FAILED: " << error.what() << '\n';
        return 1;
    }
    std::cout << "[test_fieldVectorLayout] OK\n";
    return 0;
}
