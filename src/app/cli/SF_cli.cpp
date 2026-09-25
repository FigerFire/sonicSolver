/// @file SF_cli.cpp
/// @brief 命令行解析；case 规则和系统解析均委托给 application。

#include "SF_cli.h"
#include "SF_resultCommands.h"

#include "app/application/SF_application.h"
#include "app/application/SF_inspection.h"
#include "app/application/output/SF_report.h"
#include "app/application/model/SF_configParser.h"
#include "core/config/SF_config.h"
#include "core/interfaces/SF_log.h"
#include "SF_systemPrinter.h"

#include <array>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Args = std::vector<std::string>;

void printText(const std::string& text) {
    if (SF::logFromThisProcess()) std::cout << text << std::endl;
}

void printHelp() {
    printText(R"HELP(sonicSolver - structured-grid compressible FDM solver

Usage:
  sonicSolver [CASE]                      run CASE; default CASE is .
  sonicSolver run [CASE]                 resolve, validate, then run
  sonicSolver check [CASE]               validate only; default CASE is .
  sonicSolver explain [CASE]             print the resolved equations/system
  sonicSolver doctor [CASE]              print implementation diagnostics
  sonicSolver models [CATEGORY] [OPTIONS] list available models
  sonicSolver why MODEL [--with mpi]      explain capability requirements
  sonicSolver explain-model MODEL         explain one model descriptor
  sonicSolver recipes                     list case recipes
  sonicSolver init [CASE] --recipe RECIPE  create a native case template

Commands and options:
  --initial-output [CASE]                 write initial fields; default CASE is .
  --steps N [CASE]                        limit a run to N time steps; default .
  cleanCase [CASE]                        remove CASE/result contents
  postProcess [CASE]                      open the unique result PVD
  -cleanResult [CASE]                     alias for cleanCase
  -postProcessing [CASE]                 alias for postProcess
  -h, --help                              show this help

Examples:
  sonicSolver check test/sodCase
  sonicSolver explain test/sodCase
  sonicSolver run test/sodCase
  mpirun -np 4 sonicSolver run test/IBM/IBMCase
  sonicSolver models ibm --with mpi
)HELP");
}

std::filesystem::path optionalCaseDirectory(int argc, char* argv[]) {
    if (argc > 3) {
        throw std::runtime_error("该命令最多接受一个 case 目录参数");
    }
    return argc == 3 ? std::filesystem::path(argv[2])
                     : std::filesystem::path{};
}

// 命令行 → RunRequest → application 运行入口。
//
// 保留 legacy 命令行语义（--initial-output / --steps / 默认 CASE）：
//   --initial-output CASE    初始场输出后退出
//   --steps N CASE           限制步数
//   CASE                     直接运行
// 解析在此完成；application 层只接收 RunRequest。
int runWithArgs(int argc, char* argv[]) {
    SF::Application::RunRequest request;
    request.argc = argc;
    request.argv = argv;

    if (argc < 2) {
        SF::broadcast("Usage: ", argv[0]);
        SF::broadcast("Example: ", "./sonicSolver test/sodCase");
        SF::broadcast("Initial output: ",
                      "./sonicSolver --initial-output test/sodCase");
        return 1;
    }

    int caseArgument = 1;
    if (std::string(argv[1]) == "--initial-output") {
        if (argc < 3) {
            SF::broadcast("Fatal: ",
                          "--initial-output requires a case file path.");
            return 1;
        }
        request.initialOutputOnly = true;
        caseArgument = 2;
    }

    if (std::string(argv[1]) == "--steps") {
        if (argc != 4) {
            throw std::runtime_error("Usage: sonicSolver --steps N CASE");
        }
        std::size_t consumed = 0;
        request.stepLimit = std::stoi(argv[2], &consumed);
        if (consumed != std::string(argv[2]).size()
            || request.stepLimit <= 0) {
            throw std::runtime_error("--steps requires a positive integer");
        }
        caseArgument = 3;
    }

    request.casePath = argv[caseArgument];
    return SF::Application::run(request);
}

// 命令行 - application 转换层：把逻辑参数重组为 argv 后交给 runWithArgs。
int runTranslated(const Args& arguments) {
    std::vector<std::string> translated;
    translated.reserve(arguments.size() + 1);
    translated.push_back("sonicSolver");
    for (const auto& argument : arguments) translated.push_back(argument);
    std::vector<char*> argv;
    argv.reserve(translated.size() + 1);
    for (auto& argument : translated) argv.push_back(argument.data());
    argv.push_back(nullptr);
    return runWithArgs(static_cast<int>(translated.size()), argv.data());
}

// 把 `run` 子命令的 run-control 选项重组为 legacy argv 形态，再交给
// runWithArgs 统一解析为 RunRequest。
int runCommand(int argc, char* argv[]) {
    std::string casePath;
    std::string stepLimit;
    bool initialOutput = false;
    for (int i = 2; i < argc; ++i) {
        const std::string argument = argv[i];
        if (argument == "--initial-output") {
            if (initialOutput || !stepLimit.empty())
                throw std::runtime_error("run accepts only one run-control option");
            initialOutput = true;
        } else if (argument == "--steps") {
            if (i + 1 >= argc || initialOutput || !stepLimit.empty()) {
                throw std::runtime_error("run --steps requires N CASE");
            }
            stepLimit = argv[++i];
        } else if (casePath.empty()) {
            casePath = argument;
        } else {
            throw std::runtime_error("run accepts exactly one CASE");
        }
    }
    if (casePath.empty()) casePath = ".";
    Args runArguments;
    if (initialOutput) {
        runArguments = {"--initial-output", casePath};
    } else if (!stepLimit.empty()) {
        runArguments = {"--steps", stepLimit, casePath};
    } else {
        runArguments = {casePath};
    }
    return runTranslated(runArguments);
}

const char* modelMaturity(const std::string& name) {
    if (name == "ghostCell") return "stable";
    if (name == "dfmFractionalStepPrescribed"
        || name == "velocityForcingBP") return "validated";
    return "experimental";
}

bool modelNeedsGlobalConstraint(const std::string& name) {
    return name != "ghostCell" && name != "velocityForcingBP";
}

const std::array<const char*, 10>& ibmModels() {
    static const std::array<const char*, 10> values{
        "ghostCell",
        "peskinOriginal",
        "dfmExplicitSelfPropelled",
        "dfmFractionalStepSelfPropelled",
        "dfmFractionalStepPrescribed",
        "dfmImplicitPrescribed",
        "dfmImplicitSelfPropelled",
        "dfmAugmentedLagrangian",
        "velocityForcingFTS",
        "velocityForcingBP"};
    return values;
}

void printModels(const std::string& category, bool withMPI) {
    if (!category.empty() && category != "ibm") {
        throw std::runtime_error(
            "unknown model category '" + category + "' (supported: ibm)");
    }
    printText(category.empty() ? "Available models" : "IBM models");
    if (category.empty()) {
        printText("  Flow");
        printText("    densityBased                         stable");
        printText("    pressureBased                        validated");
    }
    printText("  IBM");
    for (const char* value : ibmModels()) {
        const std::string name(value);
        const bool available = !withMPI || !modelNeedsGlobalConstraint(name);
        const std::size_t padding = name.size() < 34 ? 34 - name.size() : 1;
        printText("    " + name + std::string(padding, ' ')
                  + modelMaturity(name)
                  + (withMPI ? (available ? "  MPI=OK" : "  MPI=UNAVAILABLE")
                             : ""));
    }
}

void printWhy(const std::string& model, bool withMPI, bool detailed = false) {
    const std::string normalized = SF::FDM::normalizeToken(model);
    std::string canonical;
    for (const char* value : ibmModels()) {
        if (SF::FDM::normalizeToken(value) == normalized) {
            canonical = value;
            break;
        }
    }
    if (canonical.empty() && normalized == "ghostcell") canonical = "ghostCell";
    if (canonical.empty()) {
        throw std::runtime_error("unknown model '" + model
                                 + "'; use sonicSolver models ibm");
    }
    printText(canonical + (withMPI ? " + MPI" : ""));
    if (detailed) {
        printText("  Reference: "
                  + std::string(canonical == "ghostCell"
                                    ? "sharp-interface GhostCellIBM"
                                    : "Bhalla unified IBM formulation"));
        printText("  Strategy: "
                  + std::string(canonical == "ghostCell"
                                    ? "boundaryStencilClosure"
                                    : canonical == "velocityForcingBP"
                                        ? "dissipativePenalty"
                                        : canonical == "dfmImplicitPrescribed"
                                            || canonical == "dfmImplicitSelfPropelled"
                                            || canonical == "dfmAugmentedLagrangian"
                                            ? "monolithicKKT"
                                            : "incrementalProjection"));
        printText("  Unknowns: "
                  + std::string(canonical == "ghostCell"
                                    ? "Eulerian conservative state"
                                    : "fluid state, multiplier, optional solid DOF"));
    }
    printText("  Required capabilities:");
    printText("    ImmersedBoundaryGeometry       OK");
    if (canonical == "ghostCell") {
        printText("    ConservativeGhostState         OK");
    } else {
        printText("    PredictedConservativeState     OK");
        printText("    ImmersedConstraintProjection   OK");
        printText("    LagrangeMultiplierField        "
                  + std::string(canonical == "velocityForcingBP"
                                     ? "not-required" : "OK"));
    }
    if (withMPI) {
        const bool available = !modelNeedsGlobalConstraint(canonical);
        printText("    DistributedLinearSystem        OK");
        printText("    GlobalConstraintDof            "
                  + std::string(available ? "not-required" : "MISSING"));
        printText(available ? "  Result: AVAILABLE" :
                              "  Result: NOT AVAILABLE");
        if (!available) {
            printText("  Blocking capability: GlobalConstraintDof");
            printText("  Use a serial run or configure a supported non-multiplier method.");
        }
    } else {
        printText("  Result: AVAILABLE in the serial capability profile");
    }
}

void printRequirements(const SF::Application::CaseInspection& inspection) {
    for (const auto& requirement : inspection.system.runtime.requirements) {
        const std::string state = requirement.required
            ? (requirement.available ? "OK" : "MISSING")
            : "not-required";
        printText("  " + requirement.name + "  " + state);
    }
}

int checkCase(const std::string& path, bool detailed) {
    SF::Application::CaseInspection inspection;
    try {
        inspection = SF::Application::inspectCase(path);
    } catch (const std::exception& error) {
        printText("CONFIGURATION INVALID");
        printText("Reason: " + std::string(error.what()));
        return 1;
    }
    if (detailed) {
        printText(SF::System::describe(inspection.system)
                  + SF::IBM::renderExplain(inspection.ibmExplain));
        printText("\nCompiled plan: " + inspection.system.solvePlan.root.name);
        return 0;
    }
    printText("Checking sonicSolver configuration...");
    printText("  Flow:       equations + coupling preset  OK");
    printText("  Coupling:   "
        + std::string(SF::System::toString(inspection.system.coupling.status))
        + " (" + inspection.system.coupling.preset + ")");
    printText("  State:      " + std::string(
        SF::System::toString(inspection.system.classification.templateOrigin)) + "  OK");
    printText("  Plan:       " + inspection.system.solvePlan.root.name + "  OK");
    printText("  IBM:        "
        + std::string(inspection.ibm.enabled ? "enabled" : "disabled") + "  OK");
    printText("  Execution:");
    printRequirements(inspection);
    printText("\nCONFIGURATION VALID");
    return 0;
}

int doctorCase(const std::string& path) {
    const auto inspection = SF::Application::inspectCase(path);
    printText("[Configuration] OK");
    printText("  case: " + inspection.config.caseName);
    printText("  source: " + (inspection.config.modelDescription
        ? std::string("native case.yaml") : std::string("compatibility adapter")));
    printText("[Equation System] OK");
    printText("  compiled plan: " + inspection.system.solvePlan.root.name);
    printText("[MPI / Ownership]");
    printText("  parallel: " + std::string(
        inspection.config.parallel.enabled ? "enabled" : "disabled"));
    printRequirements(inspection);
    printText("[Output]");
    printText("  " + SF::Application::Report::formatRunControl(inspection.config));
    printText("DOCTOR COMPLETE");
    return 0;
}

void writeRecipe(const std::filesystem::path& directory,
                 const std::string& recipe) {
    const std::array<const char*, 5> recipes{
        "compressibleSingleFluid", "compressibleFixedImmersedBody",
        "compressibleMovingBody", "selfPropelledBody",
        "levelSetSurfaceTension"};
    bool known = false;
    for (const char* value : recipes) known = known || recipe == value;
    if (!known) throw std::runtime_error("unknown recipe '" + recipe + "'");
    if (std::filesystem::exists(directory)
        && !std::filesystem::is_empty(directory)) {
        throw std::runtime_error("refusing to initialize non-empty directory: "
                                 + directory.string());
    }
    std::filesystem::create_directories(directory / "solvers");
    std::filesystem::create_directories(directory / "models");
    std::filesystem::create_directories(directory / "mesh");
    std::filesystem::create_directories(directory / "fields");
    auto write = [&](const std::filesystem::path& path, const std::string& body) {
        std::ofstream output(path);
        if (!output) throw std::runtime_error("cannot write " + path.string());
        output << body;
    };
    write(directory / "case.yaml",
          "SonicFile:\n  object: case\n  type: registry\n"
          "formatVersion: 1\nname: " + directory.filename().string() + "\n");
    write(directory / "solvers/runtime.yaml",
          "SonicFile:\n  object: solver\n  type: runtime\n"
          "startTime: 0\nendTime: 1\nendStep: 100\nCFL: 0.5\n"
          "writeControl: timeStep\nwriteInterval: 10\ncreateMesh: false\n"
          "output:\n  jobName: case\n  outputDir: result\n");
    write(directory / "solvers/numerics.yaml",
          "SonicFile:\n  object: solver\n  type: numerics\n"
          "time:\n  default: forwardEuler\nterms:\n"
          "  convection: teno5Steger\n");
    write(directory / "solvers/algorithm.yaml",
          "SonicFile:\n  object: solver\n  type: algorithm\ntype: densityBase\n");
    write(directory / "solvers/solvers.yaml",
          "SonicFile:\n  object: solver\n  type: registry\n"
          "runtime:\n  type: runtime\n  file: solvers/runtime.yaml\n"
          "numerics:\n  type: numerics\n  file: solvers/numerics.yaml\n"
          "algorithm:\n  type: algorithm\n  file: solvers/algorithm.yaml\n");
    write(directory / "models/models.yaml",
          "SonicFile:\n  object: models\n  type: registry\n");
    write(directory / "mesh/mesh.yaml",
          "SonicFile:\n  object: mesh\n  type: registry\n"
          "files: [mesh/mesh.sfm]\n");
    write(directory / "fields/fields.yaml",
          "SonicFile:\n  object: fields\n  type: registry\n");
    write(directory / "fields/internalField.yaml",
          "SonicFile:\n  object: fields\n  type: internalField\n");
    write(directory / "fields/boundaries.yaml",
          "SonicFile:\n  object: fields\n  type: boundary\n");
    printText("Initialized native case template: " + directory.string());
    printText("Recipe: " + recipe);
    printText("Next: add mesh/mesh.sfm and field declarations, then run sonicSolver check "
              + directory.string());
}

} // namespace

int SF::CLI::run(int argc, char* argv[]) {
    try {
        if (argc < 2) return runTranslated({"."});
        const std::string command = argv[1];
        // 命令分发
        if (command == "-h" || command == "--help") {
            printHelp();
            return 0;
        }
        if (command == "check") {
            if (argc > 3) throw std::runtime_error(command + " accepts at most one CASE");
            return checkCase(argc == 3 ? argv[2] : ".", false);
        }
        if (command == "explain") {
            if (argc > 3) throw std::runtime_error("explain accepts at most one CASE");
            return checkCase(argc == 3 ? argv[2] : ".", true);
        }
        if (command == "doctor") {
            if (argc > 3) throw std::runtime_error("doctor accepts at most one CASE");
            return doctorCase(argc == 3 ? argv[2] : ".");
        }
        if (command == "models") {
            std::string category;
            bool withMPI = false;
            for (int i = 2; i < argc; ++i) {
                if (std::string(argv[i]) == "--with" && i + 1 < argc
                    && std::string(argv[i + 1]) == "mpi") {
                    withMPI = true;
                    ++i;
                } else if (category.empty()) category = argv[i];
                else throw std::runtime_error("models accepts one category");
            }
            printModels(category, withMPI);
            return 0;
        }
        if (command == "why") {
            if (argc < 3 || argc > 5) throw std::runtime_error(
                "why requires MODEL [--with mpi]");
            bool withMPI = false;
            for (int i = 3; i < argc; ++i) {
                if (std::string(argv[i]) == "--with" && i + 1 < argc
                    && std::string(argv[++i]) == "mpi") withMPI = true;
                else throw std::runtime_error("why supports only --with mpi");
            }
            printWhy(argv[2], withMPI);
            return 0;
        }
        if (command == "explain-model") {
            if (argc != 3) throw std::runtime_error("explain-model requires MODEL");
            printWhy(argv[2], false, true);
            return 0;
        }
        if (command == "recipes") {
            printText("Recipes:");
            printText("  compressibleSingleFluid");
            printText("  compressibleFixedImmersedBody");
            printText("  compressibleMovingBody");
            printText("  selfPropelledBody");
            printText("  levelSetSurfaceTension");
            return 0;
        }
        if (command == "init") {
            std::string casePath = ".";
            std::string recipe;
            bool caseProvided = false;
            for (int i = 2; i < argc; ++i) {
                const std::string argument = argv[i];
                if (argument == "--recipe") {
                    if (i + 1 >= argc || !recipe.empty()) {
                        throw std::runtime_error(
                            "Usage: sonicSolver init [CASE] --recipe RECIPE");
                    }
                    recipe = argv[++i];
                } else if (!caseProvided) {
                    casePath = argument;
                    caseProvided = true;
                } else {
                    throw std::runtime_error(
                        "Usage: sonicSolver init [CASE] --recipe RECIPE");
                }
            }
            if (recipe.empty()) {
                throw std::runtime_error(
                    "Usage: sonicSolver init [CASE] --recipe RECIPE");
            }
            writeRecipe(casePath, recipe);
            return 0;
        }
        if (command == "run") return runCommand(argc, argv);
        if (command == "cleanCase")
            return cleanResult(optionalCaseDirectory(argc, argv));
        if (command == "-cleanResult")
            return cleanResult(optionalCaseDirectory(argc, argv));
        if (command == "postProcess")
            return openPostProcessing(optionalCaseDirectory(argc, argv));
        if (command == "-postProcessing")
            return openPostProcessing(optionalCaseDirectory(argc, argv));
        if (command == "--initial-output" && argc == 2)
            return runTranslated({"--initial-output", "."});
        if (command == "--steps" && argc == 3)
            return runTranslated({"--steps", argv[2], "."});
        if (command == "--initial-output" || command == "--steps")
            return runWithArgs(argc, argv);

        // Allow the compact form `sonicSolver CASE`.
        if (command.front() != '-') return runWithArgs(argc, argv);
        throw std::runtime_error("unknown option '" + command
                                 + "'; use sonicSolver --help");
    } catch (const std::exception& error) {
        // MPI 非 root rank 不能只走 broadcast：若它在 collective 前失败，
        // root 会在通信中等待。每个失败进程必须直接写出自身异常，方便
        // 定位不一致的分布式控制流。
        std::cerr << "[SF] Fatal: " << error.what() << std::endl;
        SF::broadcast("Fatal: ", error.what());
        return 1;
    }
}
