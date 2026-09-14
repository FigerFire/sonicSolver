#pragma once

/// @file SF_kktDofs.h
/// @brief KKT 数学自由度身份与本地解向量槽位的独立布局。

#include "solver/linearAlgebra/SF_globalDofSystem.h"
#include "SF_immersedConstraint.h"

#include <cstdint>
#include <vector>

namespace SF::PressureBased {

/// @brief 当前网格维度真正参与 KKT 的笛卡尔速度与刚体分量。
struct KKTSubspace {
    std::vector<int> velocityComponents;
    std::vector<int> solidComponents;
};

/// @brief 从 empty 方向构造严格二维或三维 KKT 子空间。
KKTSubspace activeKKTSubspace(const FDM::ImmersedSolidEquation& solid);

/// @brief 把数学 `GlobalDofId` 与一次求解返回向量的局部槽位明确分开。
class KKTLayout {
public:
    KKTLayout(int cellCount,
              std::vector<std::int64_t> constraintEntities,
              std::vector<int> velocityComponents,
              std::vector<int> solidComponents);

    LinearAlgebra::GlobalDofId pressure(int cell) const;
    LinearAlgebra::GlobalDofId velocity(int cell, int component) const;
    LinearAlgebra::GlobalDofId constraint(int marker, int component) const;
    LinearAlgebra::GlobalDofId solid(int component) const;

    std::size_t pressureSlot(int cell) const;
    std::size_t velocitySlot(int cell, int component) const;
    std::size_t constraintSlot(int marker, int component) const;
    std::size_t solidSlot(int component) const;

    const std::vector<LinearAlgebra::GlobalDofId>& orderedDofs() const {
        return orderedDofs_;
    }

private:
    int cells_ = 0;
    std::vector<std::int64_t> constraints_;
    std::vector<int> velocityComponents_;
    std::vector<int> solidComponents_;
    std::vector<LinearAlgebra::GlobalDofId> orderedDofs_;

    void requireCell(int cell) const;
    void requireMarker(int marker) const;
    static void requireComponent(int component);
    static std::size_t componentSlot(
        const std::vector<int>& components,
        int component,
        const char* block);
};

} // namespace SF::PressureBased
