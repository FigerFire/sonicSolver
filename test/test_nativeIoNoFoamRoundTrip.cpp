/// @file test_nativeIoNoFoamRoundTrip.cpp
/// @brief §10 Gate A 回归：原生语义 YAML 必须直接产出 typed spec。
///
/// 本测试证明的是"没有中间翻译"这件事本身，而不是某个字段恰好相等：
///   1. 原生 case 目录不存在 OpenFOAM 文档布局（system/constant/0），
///      因此解码阶段不可能 round-trip 过一份 Foam 形状的中间文档；
///   2. `CaseAdapter::build` 把语义对象直接绑进 typed section，
///      field 的 type/storage/location 已经是强类型值而不是待再解析的字符串；
///   3. 完整原生路径 `ModelLoader::read` 的结果直接落在强类型 CaseConfig 上。

#include "app/application/model/SF_model.h"
#include "app/application/model/compatibility/SF_compatibility.h"
#include "infrastructure/io/case/SF_case.h"
#include "core/config/SF_runControl.h"

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>

namespace {

[[noreturn]] void fail(const std::string& message) {
    std::cerr << "native IO no-Foam round-trip test failed: " << message << '\n';
    std::exit(1);
}

void require(bool condition, const std::string& message) {
    if (!condition) fail(message);
}

} // namespace

int main() {
    const std::filesystem::path caseDir =
        std::filesystem::path(SF_TEST_SOURCE_DIR) / "test/Sod/sodCase_weno7_t0p2";
    require(std::filesystem::is_directory(caseDir),
            "Sod native case directory is missing: " + caseDir.string());

    // (1) 原生 case 里没有 OpenFOAM 文档布局可被 round-trip。
    for (const char* legacy : {"system", "constant", "0"}) {
        require(!std::filesystem::exists(caseDir / legacy),
                std::string("native case unexpectedly contains OpenFOAM layout '")
                    + legacy + "'");
    }
    require(std::filesystem::exists(caseDir / "solvers/runtime.yaml"),
            "native case does not declare solvers/runtime.yaml");

    // (2) 语义对象 -> typed section，不经过任何字典文档。
    const SF::Model::Description description =
        SF::CaseIO::read(caseDir.string());
    SF::CaseAdapter adapter(caseDir.string());
    adapter.build(description);
    const SF::NativeCaseSections& sections = adapter.sections();

    require(sections.runtime.contains("CFL"),
            "runtime semantic object was not bound into the typed runtime slot");
    require(sections.numerics.contains("terms"),
            "numerics semantic object was not bound into the typed numerics slot");
    require(sections.hasAlgorithm,
            "solver/algorithm semantic object was not recognized");
    require(sections.hasOutput && sections.hasParallel,
            "nested runtime objects were not decomposed into their own slots");

    const SF::Model::FieldDescriptor* rho = sections.field("rho");
    require(rho != nullptr, "typed field registry lost the 'rho' field");
    require(rho->type == SF::Model::ValueType::Scalar,
            "'rho' did not keep its typed scalar value type");
    require(rho->location == "point" && rho->storage == "primary",
            "'rho' did not keep its typed location/storage metadata");
    require(rho->boundaries.contains("Right")
                && rho->boundaries.contains("sideWalls"),
            "'rho' boundaries were not expanded while still in the typed model");
    const SF::Model::FieldDescriptor* velocity = sections.field("U");
    require(velocity != nullptr
                && velocity->type == SF::Model::ValueType::Vector,
            "'U' did not keep its typed vector value type");

    // (3) 完整原生路径直接产出强类型 CaseConfig。
    const SF::CaseConfig config =
        SF::Application::ModelLoader::read(caseDir.string());
    require(config.solver.numerics.timeRecipe.id()
                == SF::FDM::TimeRecipeId::ForwardEuler,
            "native time recipe was not projected onto the compiled recipe");
    require(config.time.endTime == 0.2,
            "native endTime was not projected onto RunControl");
    require(config.meshFiles == std::vector<std::string>{"mesh/mesh.sfm"},
            "native mesh file registry was not projected onto CaseConfig");
    require(config.caseName == "sodCase_weno7_t0p2",
            "native case identity was not projected onto CaseConfig");

    std::cout << "typed native IO produced CaseConfig without a Foam document\n";
    return 0;
}
