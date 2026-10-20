#pragma once

/// @file SF_algorithmDescriptor.h
/// @brief 把 IBM 用户选择解析为显式的未知量、约束和求解块描述。

#include "SF_immersedSystem.h"

namespace SF::IBM::Descriptor {

/// @brief 构造 GhostCellIBM 的边界重建描述。
FDM::ImmersedAlgorithmDescriptor ghostCell();

/// @brief 构造变分 forcing/DLM 算法的数学描述。
/// @param config 已完成 IO 解析的 IBM forcing 配置。
/// @return 不包含运行期 Field 或 MPI 状态的不可变描述。
FDM::ImmersedAlgorithmDescriptor variational(
    const FDM::IBMForcingConfig& config);

/// @brief 校验 descriptor 内部引用和配置选择的一致性。
///
/// 该校验只拒绝不完整/自相矛盾的数学系统，不会替换 enforcement 或降阶。
void validate(const FDM::ImmersedAlgorithmDescriptor& descriptor);

} // namespace SF::IBM::Descriptor
