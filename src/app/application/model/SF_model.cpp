#include "SF_model.h"
#include "compatibility/SF_compatibility.h"
#include "infrastructure/io/case/SF_case.h"
#include <filesystem>
namespace SF::Application::ModelLoader {
std::vector<EquationSystemInstanceConfig> equationInstances(
        const Model::Description& model) {
    std::vector<EquationSystemInstanceConfig> result;
    for (const auto& object : model.objects) {
        if (object.type != "equationSystem"
            && object.type != "BuiltinEquationSystem") continue;
        if (object.name.empty() || !object.parameters.contains("preset")) {
            throw std::runtime_error(
                "An equation-system instance requires a name and preset/type.");
        }
        const std::string preset =
            object.parameters.at("preset").get<std::string>();
        const std::string type = preset == "compressible"
            ? "singleFluid" : preset == "eulerianEulerian"
                ? "eulerian" : std::string{};
        if (type.empty()) {
            throw std::runtime_error(
                "Equation-system instance '"+object.name
                +"' has unsupported type '"+preset+"'.");
        }
        for (const auto& existing : result) {
            if (existing.name == object.name) {
                throw std::runtime_error(
                    "Duplicate equation-system instance name '"
                    +object.name+"'.");
            }
        }
        result.push_back({object.name,type});
    }
    return result;
}

Model::FactoryRegistry<CaseConfig>& objectRegistry() {
    static Model::FactoryRegistry<CaseConfig> registry;
    return registry;
}
CaseConfig build(const Model::Description& model) {
    const auto instances = equationInstances(model);
    auto builtins=model;
    std::vector<Model::ObjectDescriptor> extensions;
    builtins.objects.clear();
    for(const auto& object:model.objects) {
        if(objectRegistry().contains(object.type))extensions.push_back(object);
        else builtins.objects.push_back(object);
    }
    auto result=CaseAdapter(model.caseDir).build(builtins);
    for(const auto& object:extensions)objectRegistry().create(object,result);
    result.modelDescription=std::make_shared<const Model::Description>(model);
    result.composition.instances=instances;
    return result;
}
CaseConfig read(const std::string& path) {
    const std::filesystem::path p(path);
    const auto dir=std::filesystem::is_directory(p)?p:p.parent_path();
    if(std::filesystem::exists(dir/"case.yaml"))return build(CaseIO::read(path));
    throw std::runtime_error(
        "Legacy OpenFOAM case loading was removed; convert the case to case.yaml "
        "(see SF_nativeDecode.cpp / decodeNativeCase) before loading: "+path);
}
}
