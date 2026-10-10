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

#include <chrono>
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
    SF::Model::Description description =
        SF::CaseIO::read(caseDir.string());
    // The historical fixture predates mandatory native coefficients. Complete the
    // declaration in a test copy, preserving its original files and numerical values.
    for(auto& object:description.objects)if(object.type=="thermoDynamics")
        object.parameters["properties"]={{"equationOfState",{{"gamma",1.4},{"R",287.05}}},
            {"transport",{{"mu",0.},{"Pr",.72}}}};
    const auto typedCase=std::filesystem::path(SF_TEST_SOURCE_DIR)/"test/t/native-typed-authority";
    SF::CaseIO::write(description,typedCase.string());
    std::filesystem::copy_file(caseDir/"mesh/mesh.sfm",typedCase/"mesh/mesh.sfm",std::filesystem::copy_options::overwrite_existing);
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
        SF::Application::ModelLoader::read(typedCase.string());
    require(config.solver.numerics.timeRecipe.id()
                == SF::FDM::TimeRecipeId::ForwardEuler,
            "native time recipe was not projected onto the compiled recipe");
    require(config.time.endTime == 0.2,
            "native endTime was not projected onto RunControl");
    require(config.meshFiles == std::vector<std::string>{"mesh/mesh.sfm"},
            "native mesh file registry was not projected onto CaseConfig");
    require(config.caseName == "sodCase_weno7_t0p2",
            "native case identity was not projected onto CaseConfig");

    require(config.composition.stateDeclared
        && config.composition.solutionVariables==std::vector<std::string>{"rho","rhoU","rhoE"}
        && config.composition.algorithm=="Explicit",
        "Native WHAT/STATE/HOW declarations were replaced by compatibility inference");
    auto withoutState=description;
    withoutState.objects.erase(std::remove_if(withoutState.objects.begin(),withoutState.objects.end(),
        [](const auto& object) { return object.type=="stateRegistry"; }),withoutState.objects.end());
    bool missingStateRejected=false;
    try { (void)SF::CaseAdapter(caseDir.string()).build(withoutState); }
    catch (const std::runtime_error& error) {
        missingStateRejected=std::string(error.what()).find("explicit solution STATE")!=std::string::npos;
    }
    require(missingStateRejected,"A native EOS or legacy label silently supplied missing STATE");
    for(const auto* wall:{"noSlip","slip"}) {
        auto model=description;
        auto u=std::find_if(model.fields.begin(),model.fields.end(),[](const auto& f){return f.name=="U";});
        u->boundaries["sideWalls"]={{"type",wall}};
        SF::Model::FieldDescriptor temperature;temperature.name="T";temperature.initial=300.;
        temperature.boundaries["sideWalls"]={{"type","fixedTemperature"},{"value",310.}};
        model.fields.push_back(temperature);
        const auto resolved=SF::Application::ModelLoader::build(model);
        const auto& velocity=resolved.solver.boundaries.velocity;
        const auto bc=std::find_if(velocity.begin(),velocity.end(),[](const auto& x){return x.name=="sideWalls";});
        require(bc!=velocity.end() && bc->type==(std::string(wall)=="noSlip"?SF::FIXED_VALUE:SF::SYMMETRY)
            && bc->value.x==0. && bc->value.y==0. && bc->value.z==0.,"Wall preset did not expand to its neutral velocity constraint");
        const auto& thermal=resolved.solver.boundaries.thermal;
        require(thermal.size()==1 && thermal.front().type==SF::ThermalBCType::FixedTemperature
            && thermal.front().value==310.,"Wall velocity silently replaced the independently selected thermal law");
        require(resolved.solver.boundaries.pressure.size()==config.solver.boundaries.pressure.size()
            && resolved.solver.numerics.timeRecipe.id()==config.solver.numerics.timeRecipe.id(),"Wall preset changed pressure/time composition");
        auto invalid=model;
        auto rho=std::find_if(invalid.fields.begin(),invalid.fields.end(),[](const auto& f){return f.name=="rho";});
        rho->boundaries["sideWalls"]={{"type",wall}};
        bool rejected=false;try{(void)SF::Application::ModelLoader::build(invalid);}catch(const std::runtime_error& e){rejected=std::string(e.what()).find("velocity patches")!=std::string::npos;}
        require(rejected,"Scalar field accepted a velocity wall law");
        invalid=model;invalid.fields.back().boundaries["sideWalls"]["type"]=wall;
        rejected=false;try{(void)SF::Application::ModelLoader::build(invalid);}catch(const std::runtime_error& e){rejected=std::string(e.what()).find("thermal boundary independently")!=std::string::npos;}
        require(rejected,"Thermal field accepted a velocity wall law");
        invalid=model;auto moving=std::find_if(invalid.fields.begin(),invalid.fields.end(),[](const auto& f){return f.name=="U";});
        moving->boundaries["sideWalls"]["value"]=SF::Model::Parameters::array({1.,0.,0.});
        rejected=false;try{(void)SF::Application::ModelLoader::build(invalid);}catch(const std::runtime_error& e){rejected=std::string(e.what()).find("stationary")!=std::string::npos;}
        require(rejected,"Stationary wall ignored a conflicting velocity");
    }
    auto unknown=description;auto unknownU=std::find_if(unknown.fields.begin(),unknown.fields.end(),[](const auto& f){return f.name=="U";});
    unknownU->boundaries["sideWalls"]={{"type","slipp"}};
    bool unknownRejected=false;try{(void)SF::Application::ModelLoader::build(unknown);}catch(const std::runtime_error& e){unknownRejected=std::string(e.what()).find("Unknown boundary law")!=std::string::npos;}
    require(unknownRejected,"Unknown boundary law silently fell back");
    const auto roundTrip=std::filesystem::path(SF_TEST_SOURCE_DIR)/"test/t"/
        ("sonic-explicit-state-io-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    SF::CaseIO::write(description,roundTrip.string());
    const auto reread=SF::CaseIO::read(roundTrip.string());
    for (const auto* type:{"stateRegistry","equationRegistry","algorithmRegistry"}) {
        const auto find=[&](const auto& model) {
            return std::find_if(model.objects.begin(),model.objects.end(),
                [&](const auto& object) { return object.type==type; });
        };
        const auto before=find(description),after=find(reread);
        require(before!=description.objects.end() && after!=reread.objects.end()
            && before->parameters==after->parameters,
            std::string("Native registry round trip changed ")+type);
    }
    std::filesystem::remove_all(roundTrip);

    std::cout << "typed native IO produced CaseConfig without a Foam document\n";
    return 0;
}
