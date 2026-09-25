#pragma once

/// @file SF_globalEntityId.h
/// @brief IBM 与分布式装配共享的强类型全局实体身份。
/// @details 这些值是数学/拓扑身份，不是 patch-local index、MPI rank 或线性后端行号。

#include <cstdint>
#include <stdexcept>

namespace SF {

/// @brief 稳定的 Lagrangian 表面 marker 身份。
class GlobalMarkerId {
public:
    GlobalMarkerId() = default;

    static GlobalMarkerId fromSurfacePrimitive(std::int64_t primitive) {
        if (primitive < 0) {
            throw std::runtime_error("GlobalMarkerId requires a non-negative surface primitive.");
        }
        return GlobalMarkerId(primitive);
    }

    std::int64_t value() const { return value_; }
    bool valid() const { return value_ >= 0; }

    friend bool operator==(GlobalMarkerId a, GlobalMarkerId b) {
        return a.value_ == b.value_;
    }
    friend bool operator!=(GlobalMarkerId a, GlobalMarkerId b) {
        return !(a == b);
    }

private:
    explicit GlobalMarkerId(std::int64_t value) : value_(value) {}
    std::int64_t value_ = -1;
};

/// @brief 表面乘子/无滑移约束自由度的身份；与 marker 和 backend row 分离。
class GlobalConstraintDofId {
public:
    GlobalConstraintDofId() = default;

    static GlobalConstraintDofId fromMarker(GlobalMarkerId marker) {
        if (!marker.valid()) {
            throw std::runtime_error("GlobalConstraintDofId requires a valid marker.");
        }
        return GlobalConstraintDofId(marker.value());
    }

    std::int64_t value() const { return value_; }
    bool valid() const { return value_ >= 0; }

    friend bool operator==(GlobalConstraintDofId a, GlobalConstraintDofId b) {
        return a.value_ == b.value_;
    }
    friend bool operator!=(GlobalConstraintDofId a, GlobalConstraintDofId b) {
        return !(a == b);
    }

private:
    explicit GlobalConstraintDofId(std::int64_t value) : value_(value) {}
    std::int64_t value_ = -1;
};

} // namespace SF
