/// @file SF_hostCapabilities.cpp
/// @brief host backend 能力检测。
///
/// Data flow:
///   linked backends (infrastructure / linear algebra)
///       -> BuildCapabilities
///       -> BuildRequest
///       -> RuntimeRequirement
///
/// 这里只汇总“本 binary 实际链接了什么”。case 需要什么由 equation/system
/// 编译结果回答；两者在 RuntimeRequirement 里相遇。

#include "SF_hostCapabilities.h"

// 两个 backend 通过互斥实现回答自身是否真的存在：
// - Parallel::mpiCompiled()   : SF_mpi 真实后端 vs serial stub
// - LinearAlgebra::hypreBackendLinked() : 真实 HYPRE ParCSR vs 占位实现
#include "SF_parallelContext.h"
#include "solver/linearAlgebra/hypre/SF_hypre.h"

namespace SF::Application {

System::BuildCapabilities detectBuildCapabilities() {
    System::BuildCapabilities capabilities;
    capabilities.mpi = Parallel::mpiCompiled();
    capabilities.hypre = LinearAlgebra::hypreBackendLinked();
    // 分布式线性系统 = GlobalDofId 到 ParCSR row 的映射；它需要同时存在
    // 真实通信后端和 HYPRE 后端。
    capabilities.distributedLinearSystem =
        capabilities.mpi && capabilities.hypre;
    // canonical ConstraintGlobalDof ownership 依赖真实 owner->COPY 通信。
    // serial stub 的 copyCanonical 是 no-op，不能构成该能力。
    capabilities.canonicalConstraintDof = capabilities.mpi;
    return capabilities;
}

} // namespace SF::Application
