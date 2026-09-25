#pragma once

/// @file SF_environment.h
/// @brief execution environment 的构建与 Phase-2 runtime validation。
///
/// Program（SF_resolvedSimulationSystem）描述“解什么、怎么排”；
/// Environment 描述“在什么 mesh / storage / MPI 拓扑上解”。
/// 数学编译保持 execution-topology 无关，这里的对象只负责把已经编译好的
/// program 落到具体的存储与并行布局上。

#include "SF_caseConfig.h"
#include "SF_compositeIBM.h"
#include "SF_field.h"
#include "SF_IBM.h"
#include "SF_MultiBlockMesh.h"
#include "SF_parallelContext.h"
#include "solver/system/SF_runtimeRequirements.h"

#include <memory>

namespace SF {

class ResultWriter;

namespace Application {

struct CaseInspection;

/// @brief storage / topology 布局，由 DomainBuilder 根据 mesh + MPI 决定。
///
/// 它回答“数据在哪里、怎么同步”，不回答“解什么方程”。
enum class DomainKind {
    /// 串行 Field，或 MPI 每 rank 单 structured patch 的 Field。
    Single,
    /// MPI 多 patch 复合网格，由 executeMulti 适配。
    DistributedMultiPatch,
};

/// @brief 一次运行的执行环境。
///
/// 所有权规则：
/// - `parallel` 是 move-only 的 ParallelContext，由 build() 创建。
/// - `field` / `parallelMesh` / `compositeIBM` / `ibm` 是运行时资源，
///   其生命周期覆盖整个 run()，因此 activeField 指向其中的
///   Field 在本次运行内保持稳定。
/// - `writer` 不拥有 ResultWriter（由 run() 栈持有）。
struct ExecutionEnvironment {
    std::unique_ptr<Parallel::ParallelContext> parallel;

    /// 串行 Field 与 MPI 单 patch 统一由 activeField 指向；
    /// 多 patch 由 parallelMesh 承载。
    Field field;
    Field* activeField = nullptr;
    MultiBlockMesh parallelMesh;

    IBM::CompositeIB compositeIBM;
    IBM::IB ibm;

    ResultWriter* writer = nullptr;
    DomainKind domain = DomainKind::Single;

    bool ibmEnabled = false;
    int localBlockId = 0;

    /// IBM setup / method 是 distributed precondition：本地结果先记录，
    /// 由 validate() 统一做 allRanksAgree 归约。
    bool localIBMSetupOk = true;
    bool ibmUsesForcing = false;
};

namespace Environment {

/// @brief Phase-1 通过后构建执行环境。
///
/// 只负责 mesh / MPI / storage / IBM runtime 资源与本地（非 collective）
/// 检查；失败的 Fatal 信息在此输出。
/// @return false 时调用方直接返回非零退出码。
bool build(ExecutionEnvironment& env,
           const CaseConfig& caseConfig,
           const CaseInspection& inspection,
           const System::RuntimeRequirements& requirements,
           int& argc, char**& argv,
           ResultWriter& writer);

/// @brief Phase-2 校验：只做依赖实际 runtime topology 的 collective 检查。
///
/// 例如 IBM setup / method 是否在所有 rank 上一致——这些在没有 MPI 环境时
/// 无法判断，因此与 application 的 validateProgram 分开。
/// @return false 时已输出 Fatal。
bool validate(ExecutionEnvironment& env);

/// @brief createMesh 模式：只生成网格，不求解。
/// @return 进程退出码（0 / -1）。
int buildMeshOnly(const CaseConfig& caseConfig,
                  int& argc, char**& argv,
                  ResultWriter& writer);

} // namespace Environment
} // namespace Application
} // namespace SF
