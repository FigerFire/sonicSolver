/// @file SF_run.cpp
/// @brief 应用层 composition root：parse → read → compile → validate → build → validate → execute。

/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.03.23-----------*/

#include "SF_run.h"
#include "SF_environment.h"
#include "SF_preflight.h"
#include "app/application/execution/SF_runners.h"
#include "SF_resultWriter.h"
#include "app/application/model/SF_model.h"
#include "app/application/model/SF_output.h"
#include "app/application/system/SF_inspection.h"
#include "SF_systemPrinter.h"
#include "core/interfaces/SF_log.h"

#include <algorithm>
#include <exception>
#include <optional>
#include <string>

namespace SF::Application {
namespace {

/// @brief Run::execute —— 只跑 CompiledSolvePlan；storage 布局由 env.domain 决定。
///
/// single / multi 只存在于 execution environment 的拓扑里，不再暴露成
/// application runner 的两个独立入口名。
int execute(const System::ResolvedSimulationSystem& system,
            const System::CompiledSolvePlan& plan,
            const CaseConfig& caseConfig,
            ExecutionEnvironment& env,
            const Preflight::RunRequest& request) {
    if (env.domain == DomainKind::DistributedMultiPatch) {
        return Runners::runMultiPatch(
            env.parallelMesh, *env.writer, *env.parallel, env.compositeIBM,
            caseConfig.solver, system, plan, caseConfig, env.ibmEnabled,
            request.initialOutputOnly);
    }
    return Runners::runSingleField(
        *env.activeField, *env.writer, *env.parallel, env.ibm,
        caseConfig.solver, system, plan, caseConfig, env.ibmEnabled,
        request.initialOutputOnly, env.localBlockId);
}

} // namespace

// Application 层真正执行一个算例的主入口
int runConfiguredCase(int argc, char* argv[]) {
    const Preflight::RunRequest request =
        Preflight::parseCommandLine(argc, argv);
    if (request.usageRequested) {
        return Preflight::reportUsage(request, argv[0]);
    }
    if (!request.ok) {
        SF::broadcast("Fatal: ", request.error);
        return request.exitCode;
    }

    CaseConfig caseConfig = ModelLoader::read(request.casePath);
    if (request.stepLimit > 0) {
        caseConfig.time.endStep = caseConfig.time.endStep > 0
            ? std::min(request.stepLimit, caseConfig.time.endStep)
            : request.stepLimit;
    }

    ResultWriter resultWriter(caseConfig.output);
    ModelLoader::configureOutput(resultWriter, caseConfig);

    if (caseConfig.createMesh) {
        return Environment::buildMeshOnly(caseConfig, argc, argv,
                                          resultWriter);
    }

    // 编译 mathematical system：完全不依赖 MPI rank / partition / single-multi。
    // 编译失败（如方程/变换/能力组合不完整）属于 workflow 级兼容性问题，
    // 与旧实现一致：Fatal 后返回 -1，而不是把异常抛给 CLI。
    std::optional<CaseInspection> inspection;
    try {
        inspection.emplace(inspectCase(caseConfig));
    } catch (const std::exception& error) {
        SF::broadcast("Fatal workflow compatibility: ", error.what());
        return -1;
    }

    SF::broadcast("Compiled plan     : ",
                  inspection->system.solvePlan.root.name);
    SF::broadcast("", System::describe(inspection->system));

    if (!Preflight::validateProgram(inspection->system)) {
        return -1;
    }

    ExecutionEnvironment env;
    if (!Environment::build(env, caseConfig, *inspection, argc, argv,
                            resultWriter)) {
        return -1;
    }
    if (!Environment::validate(env)) {
        return -1;
    }

    return execute(inspection->system, inspection->system.solvePlan,
                   caseConfig, env, request);
}

} // namespace SF::Application
