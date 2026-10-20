#pragma once

/// @file SF_globalDofResidual.h
/// @brief GlobalDof 残差的 owner 汇总与副本广播。

#include "SF_communicationPlan.h"
#include "core/residual/SF_residual.h"

#include <vector>

namespace SF {
struct MeshBlockField;
class Field;

namespace Parallel {
namespace Backend { class CommunicationBackend; }

/// @brief 把各结构 patch 的积分残差汇总为唯一 GlobalDof 残差。
class GlobalDofResidualAssembler {
public:
    /// @brief 执行 `contributors -> owner sum -> owner broadcast`。
    /// @param plan 预计算的 GlobalDof ownership 与对偶体积片段。
    /// @param blocks 全局 patch 容器；每个 rank 只读取自己拥有的 patch。
    /// @param backend 可为空；为空时执行确定性的单 rank 装配。
    static void assemble(
        const MeshCommunication::HaloExchangePlan& plan,
        std::vector<MeshBlockField>& blocks,
        const std::vector<Residual*>& residuals,
        const Backend::CommunicationBackend* backend);

    /// @brief 一块/每 rank 布局的同一装配语义。
    static void assemble(
        const MeshCommunication::HaloExchangePlan& plan,
        Field& localField,
        Residual& residual,
        int localBlockId,
        const Backend::CommunicationBackend* backend);
};

} // namespace Parallel
} // namespace SF
