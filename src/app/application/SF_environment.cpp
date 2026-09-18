/// @file SF_environment.cpp
/// @brief execution environment 的构建与 runtime validation。

#include "SF_environment.h"

#include "SF_boundaryGeometry.h"
#include "SF_mesh.h"
#include "app/application/output/SF_fields.h"
#include "SF_resultWriter.h"
#include "app/application/model/SF_runtimeConfig.h"
#include "app/application/system/SF_inspection.h"

#include "core/interfaces/SF_log.h"

#include <string>
#include <vector>

namespace SF::Application::Environment {
namespace {

/// @brief DomainBuilder：由 mesh topology + MPI 决定 storage 布局。
///
/// 数学 program 不参与。每个 rank 都会独立持有一份 partition metadata；
/// 路由选择必须先归约成一个全局决定，否则某些 rank 进入 single-field 的
/// collective、另一些 rank 提前走 multi-patch / 结束路径，会形成无法诊断的
/// MPI 死锁。
bool resolveSingleFieldLayout(const MultiBlockMesh& mesh,
                              Parallel::ParallelContext& parallel) {
    const bool localMultiFieldLayout =
        Output::needsMultiFieldMPI(mesh, parallel.size());
    return parallel.allRanksAgree(!localMultiFieldLayout);
}

} // namespace

bool build(ExecutionEnvironment& env,
           const CaseConfig& caseConfig,
           const CaseInspection& inspection,
           int& argc, char**& argv,
           ResultWriter& writer) {
    const FDM::SolverConfig& solverConfig = caseConfig.solver;
    const MeshRuntimeConfig& meshConfig = inspection.mesh;
    const IBM::IBMRuntimeConfig& ibmConfig = inspection.ibm;

    // Pressure-based HYPRE backend 直接使用 MPI_COMM_WORLD 并调用
    // HYPRE_Initialize；即使 parallel.enabled=false，单 rank 也必须先
    // MPI_Init 才能让 HYPRE 工作。因此 distributed context 的请求条件是
    // “并行已启用 或 线性代数 provider 需要 MPI-initialized context”，
    // 而不是旧代码里 solver == PressureBased 这一身份判断。
    const bool needsDistributedContext =
        caseConfig.parallel.enabled
        || solverConfig.numerics.solver
               == FDM::SolverAlgorithm::PressureBased;

    env.parallel = std::make_unique<Parallel::ParallelContext>(
        argc, argv, needsDistributedContext);
    writer.setParallelCoordinator(&env.parallel->coordinator());
    env.writer = &writer;

    if (caseConfig.parallel.enabled && !env.parallel->available()) {
        SF::broadcast("Fatal: ", "This build was not linked with MPI.");
        return false;
    }

    if (caseConfig.parallel.enabled && env.parallel->active()
        && caseConfig.parallel.processCount != env.parallel->size()) {
        if (env.parallel->isRoot()) {
            SF::broadcast("Fatal: ",
                          "runtime MPI ranks="
                              + std::to_string(env.parallel->size())
                              + " must equal [parallel].split product="
                              + std::to_string(caseConfig.parallel.processCount));
        }
        return false;
    }

    // variational forcing IBM 的选型信息在配置期即可确定，与 runtime 无关。
    if (ibmConfig.enabled
        && ibmConfig.method == FDM::IBMMethod::VariationalForcing) {
        SF::broadcast(
            "IBM selection     : ",
            std::string("algorithm=")
                + FDM::toString(ibmConfig.forcing.algorithm)
                + ", support="
                + FDM::toString(ibmConfig.forcing.constraintSupport)
                + ", representation="
                + FDM::toString(ibmConfig.forcing.representation)
                + ", enforcement="
                + FDM::toString(ibmConfig.forcing.enforcement)
                + ", solid="
                + FDM::toString(ibmConfig.forcing.solidModel));
    }

    env.ibmEnabled = ibmConfig.enabled;
    env.compositeIBM.configure(ibmConfig);

    if (caseConfig.parallel.enabled && env.parallel->active()) {
        MultiBlockMesh::CompositePreprocessor preprocessor;
        // 只有 sharp-interface Ghost IBM 需要在 MPI 切分前写入 source-zone
        // SDF/ghost 分类。Variational forcing 保持 Eulerian 流体点为
        // canonical unknown，并在每个 rank 的 IB facade 内构造 owner-local
        // J edge。
        if (ibmConfig.enabled
            && ibmConfig.method == FDM::IBMMethod::Ghost) {
            preprocessor = [&](std::vector<RawMeshBlock>& zones) {
                return env.compositeIBM.preprocessSourceZones(
                    zones, caseConfig.caseDir,
                    caseConfig.ibmGeometryFiles);
            };
        }
        if (!Mesh::setupMultiBlockMesh(
                env.parallelMesh, caseConfig, meshConfig, preprocessor)) {
            return false;
        }
        if ((int)env.parallelMesh.partitionCount() != env.parallel->size()) {
            if (env.parallel->isRoot()) {
                SF::broadcast("Fatal: ",
                              "decompose produced "
                                  + std::to_string(env.parallelMesh.partitionCount())
                                  + " partitions for "
                                  + std::to_string(env.parallel->size())
                                  + " MPI ranks.");
            }
            return false;
        }

        writer.configurePieces(
            Output::buildMultiBlockVTKPieces(env.parallelMesh));

        // MultiPatch runner 当前是 Ghost boundary 的复合 patch 适配器。对于
        // variational forcing，单 source-zone 的每 rank 单 patch 直接进入统一
        // Field/ExecutionRuntime 路径，才能让 surface J/J^T/lambda 使用新
        // ConstraintGlobalDof ownership；多 patch forcing 仍明确拒绝而非伪装。
        const bool singleFieldLayout =
            resolveSingleFieldLayout(env.parallelMesh, *env.parallel);
        const bool forcingIBM = ibmConfig.enabled
            && ibmConfig.method == FDM::IBMMethod::VariationalForcing;
        if (forcingIBM && !singleFieldLayout) {
            if (env.parallel->isRoot()) {
                SF::broadcast(
                    "Fatal: ",
                    "variational IBM MPI requires one local structured "
                    "Field per rank; multi-patch constraint routing is "
                    "not implemented.");
            }
            return false;
        }

        const bool useMultiFieldMPI =
            (ibmConfig.enabled && ibmConfig.method == FDM::IBMMethod::Ghost)
            || !singleFieldLayout;
        if (useMultiFieldMPI) {
            env.domain = DomainKind::DistributedMultiPatch;
            // runMultiPatch 自行完成 empty-dimensions 与 setupLocalPatches；
            // 此处不执行 single-field 的 IBM collective 检查。
            return true;
        }

        if ((int)env.parallelMesh.size() != env.parallel->size()) {
            if (env.parallel->isRoot()) {
                SF::broadcast("Fatal: ",
                              "MPI solve currently requires one mesh block per rank. blocks="
                                  + std::to_string(env.parallelMesh.size())
                                  + ", ranks="
                                  + std::to_string(env.parallel->size()));
            }
            return false;
        }

        env.localBlockId = env.parallel->rank();
        env.activeField =
            &env.parallelMesh.block((size_t)env.localBlockId).field;
        const bool canonicalInterfaceFlux =
            solverConfig.numerics.interfaceFlux
            == FDM::InterfaceFluxPolicy::SharedInterfaceFlux;
        env.parallel->configureSingle(
            &env.parallelMesh.haloExchangePlan(), env.localBlockId,
            &env.parallelMesh.blocks(), canonicalInterfaceFlux);
        // 保守量尚未由 EquationSet 从用户初值闭合。初始 halo 必须等
        // StateBundle 注册后由 ExecutionRuntime 的显式 ReadHalo contract
        // 触发；此处交换会把零初始化 payload 错当作 canonical 流场状态。
    } else {
        if (!Mesh::setupComplexMesh(env.field, caseConfig, meshConfig)) {
            return false;
        }
        env.activeField = &env.field;
        env.domain = DomainKind::Single;
    }

    // 边界空维度配置（single-field 路径；multi-patch 由 runner 自己配置）。
    // 先配置空维度，避免在 Field 内部构造无效索引。
    Boundary::Geometry::configureEmptyDimensions(
        *env.activeField,
        solverConfig.boundaries.density,
        solverConfig.boundaries.velocity,
        solverConfig.boundaries.energyFromPressure,
        solverConfig.boundaries.thermal,
        solverConfig.turbulence.scalars.kBoundary,
        solverConfig.turbulence.scalars.epsilonBoundary,
        solverConfig.turbulence.scalars.omegaBoundary);

    // IBM 本地 setup；collective 校验在 validate() 统一处理。
    if (env.ibmEnabled) {
        SF::broadcast("IBM enabled: ", caseConfig.ibmGeometryFiles.size());
        env.localIBMSetupOk = env.ibm.setup(
            *env.activeField, caseConfig.caseDir,
            caseConfig.ibmGeometryFiles, ibmConfig);
        env.ibmUsesForcing = env.ibm.usesForcing();
    }

    return true;
}

bool validate(ExecutionEnvironment& env) {
    // multi-patch runner 在 runMultiPatch 内部完成 setupLocalPatches 与错误
    // 处理，此处只校验 single-field 路径的 IBM distributed precondition。
    if (env.domain == DomainKind::DistributedMultiPatch) {
        return true;
    }
    if (!env.ibmEnabled) {
        return true;
    }

    // IBM geometry/method initialization is a distributed precondition.
    // No rank may silently continue without IBM while another rank enters
    // a surface J/J^T collective; that would mismatch MPI collective order.
    if (!env.parallel->allRanksAgree(env.localIBMSetupOk)) {
        if (env.parallel->isRoot()) {
            SF::broadcast("Fatal: ",
                          "IBM setup failed on at least one MPI rank.");
        }
        return false;
    }
    const bool allForcing =
        env.parallel->allRanksAgree(env.ibmUsesForcing);
    const bool allNonForcing =
        env.parallel->allRanksAgree(!env.ibmUsesForcing);
    if (!allForcing && !allNonForcing) {
        if (env.parallel->isRoot()) {
            SF::broadcast("Fatal: ",
                          "IBM method selection differs between MPI ranks.");
        }
        return false;
    }
    return true;
}

int buildMeshOnly(const CaseConfig& caseConfig,
                  int& argc, char**& argv,
                  ResultWriter& writer) {
    const MeshRuntimeConfig meshConfig =
        Runtime::makeMeshRuntimeConfig(caseConfig, caseConfig.solver);

    Parallel::ParallelContext parallel(argc, argv, false);
    writer.setParallelCoordinator(&parallel.coordinator());
    if (caseConfig.parallel.enabled && !parallel.available()) {
        SF::broadcast("Fatal: ", "This build was not linked with MPI.");
        return -1;
    }

    bool ok = true;
    if (!parallel.active() || parallel.isRoot()) {
        ok = Mesh::createMesh(caseConfig, meshConfig);
    }
    return parallel.allRanksAgree(ok) ? 0 : -1;
}

} // namespace SF::Application::Environment
