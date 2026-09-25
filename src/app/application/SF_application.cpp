/*--------------Sonic Fluid-------------------*/
/*---------copyright by Li Pengfei------------*/
/*----------code by LPF, 2026.07.30-----------*/

// Application::run
//  │
//  ├─ 1. ModelLoader::read
//  │      ↓
//  │    CaseConfig
//  │
//  ├─ 2. inspectCase
//  │      ↓
//  │    CaseInspection
//  │      └─ ResolvedSimulationSystem
//  │
//  ├─ 3. Environment::build
//  │      ↓
//  │    ExecutionEnvironment
//  │
//  └─ 4. Execution::execute
//         ↓
//       真正开始计算

#include "SF_application.h"
#include "SF_environment.h"
#include "app/application/execution/SF_execution.h"
#include "SF_resultWriter.h"
#include "app/application/model/SF_model.h"
#include "app/application/model/SF_output.h"
#include "app/application/SF_inspection.h"
#include "SF_systemPrinter.h"
#include "core/interfaces/SF_log.h"

#include <algorithm>
#include <exception>
#include <optional>
#include <string>

namespace SF::Application {
namespace {

/// @brief Phase-1 fail-fast：只验证数学/算法 program 的可执行能力。
/// @return true 表示可进入环境构建；false 时已输出 Fatal 原因。
bool validateProgram(const System::ResolvedSimulationSystem& system) {
    if (system.runtime.report.status == System::RuntimeStatus::Unsupported) {
        SF::broadcast(
            "Fatal execution capability: ",
            system.runtime.report.reason.empty()
                ? "resolved mathematical system has no supported execution path"
                : system.runtime.report.reason);
        return false;
    }
    return true;
}

} // namespace

// Application 层真正执行一个算例的主入口
int run(const RunRequest& request) {
    CaseConfig caseConfig = ModelLoader::read(request.casePath);
    if (request.stepLimit > 0) {
        caseConfig.time.endStep = caseConfig.time.endStep > 0
            ? std::min(request.stepLimit, caseConfig.time.endStep)
            : request.stepLimit;
    }

    ResultWriter resultWriter(caseConfig.output);
    ModelLoader::configureOutput(resultWriter, caseConfig);

    int runtimeArgc = request.argc;
    char** runtimeArgv = request.argv;
    if (caseConfig.createMesh) {
        return Environment::buildMeshOnly(caseConfig, runtimeArgc, runtimeArgv,
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
    SF::broadcast("", System::describe(inspection->system)
                    + IBM::renderExplain(inspection->ibmExplain));

    if (!validateProgram(inspection->system)) {
        return -1;
    }

    ExecutionEnvironment env;
    if (!Environment::build(env, caseConfig, *inspection,
                            inspection->system.runtime,
                            runtimeArgc, runtimeArgv,
                            resultWriter)) {
        return -1;
    }
    if (!Environment::validate(env)) {
        return -1;
    }

    return Execution::execute(inspection->system, inspection->system.solvePlan,
                              caseConfig, env, request);
}

} // namespace SF::Application
