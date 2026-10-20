#pragma once

/// @file SF_buildCapabilities.h
/// @brief 本 binary 真实提供的能力集合。
///
/// 该 value object 描述“这个可执行文件里链接了哪些 backend”，而不是 case
/// 需要什么。它是编译/运行边界上唯一的 backend 事实来源：equation/system
/// 编译器把 case 需求与这份 host 事实相比，得到 RuntimeRequirement。
///
/// 禁止在数学系统构造处写 `constexpr bool ... = true`：那是把构建配置伪装成
/// 数学事实，会让缺少 backend 的 binary 在进入 timestep 后才崩溃。

#include <string>

namespace SF::System {

/// @brief 由真实链接/注册结果得出的 host 能力。
struct BuildCapabilities {
    /// @brief 本 binary 是否链接了真实 MPI 通信后端。
    bool mpi = false;
    /// @brief 本 binary 是否链接了 HYPRE ParCSR 后端。
    bool hypre = false;
    /// @brief 本 binary 是否能把 GlobalDofId 映射到分布式线性系统 row。
    bool distributedLinearSystem = false;
    /// @brief 本 binary 是否提供 canonical ConstraintGlobalDof ownership/COPY。
    bool canonicalConstraintDof = false;

    std::string describe() const;
};

inline std::string BuildCapabilities::describe() const {
    std::string result;
    const auto add = [&result](bool present, const char* name) {
        if (!present) return;
        if (!result.empty()) result += ",";
        result += name;
    };
    add(mpi, "mpi");
    add(hypre, "hypre");
    add(distributedLinearSystem, "distributedLinearSystem");
    add(canonicalConstraintDof, "canonicalConstraintDof");
    return result.empty() ? "serial-only" : result;
}

} // namespace SF::System
