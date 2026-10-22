#include "SF_model.h"
#include "compatibility/SF_compatibility.h"
#include "infrastructure/io/case/SF_case.h"
#include <filesystem>
#include <cctype>
namespace SF::Application::ModelLoader {
namespace {
using P=Model::Parameters;
Model::Description expandPresets(const Model::Description& input) {
    auto result=input;
    std::map<std::string,Model::ObjectDescriptor> defaults;
    P numerical=P::object();
    std::map<std::string,std::string> origins;
    const auto role=[](std::string name) {
        for (auto& c:name)c=static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        if(name=="builtinequationsystem")name="equationsystem";
        return name;
    };
    const auto explicitlySet=[&](const std::string& path) {
        const std::string numericalPrefix="solver/numerics";
        if(path.rfind(numericalPrefix+"/",0)==0)
            return input.numerics.contains(P::json_pointer(path.substr(numericalPrefix.size())));
        for(const auto& object:input.objects) {
            const auto prefix="object/"+role(object.type);
            if(object.type!="compositionPreset" && path.rfind(prefix+"/",0)==0
                && object.parameters.contains(P::json_pointer(path.substr(prefix.size()))))return true;
        }
        return false;
    };
    const auto record=[&](const std::string& path,const std::string& source,bool explicitValue) {
        const auto old=origins.find(path);
        const auto previous=old==origins.end()?std::string():old->second;
        origins[path]=(explicitValue?"explicit ":"default ")+source;
        if(explicitValue && !previous.empty())origins[path]+="; overrides "+previous;
        else if(!explicitValue && !previous.empty() && previous!=origins[path])origins[path]=previous+"; also default "+source;
    };
    const auto merge=[&](const auto& self,P& target,const P& value,const std::string& path,
                         const std::string& source,bool explicitValue)->void {
        if (value.is_object()) {
            if (!target.is_object()) {
                if (!target.is_null() && !target.empty()) throw std::runtime_error("Conflicting setting shape at "+path+": "+origins[path]+" and "+source);
                target=P::object();
            }
            if(value.empty())record(path,source,explicitValue);
            for (const auto& item:value.items()) self(self,target[item.key()],item.value(),path+"/"+item.key(),source,explicitValue);
            return;
        }
        if (!target.is_null() && target!=value && !explicitValue && !explicitlySet(path))
            throw std::runtime_error("Ambiguous preset default at "+path+": "+origins[path]+" and "+source);
        target=value;record(path,source,explicitValue);
    };
    std::vector<Model::ObjectDescriptor> explicitObjects;
    for (const auto& object:input.objects) {
        if (object.type!="compositionPreset") {explicitObjects.push_back(object);continue;}
        Model::Schema{{{"preset","string",true}},false}.validate(object.parameters,object.source);
        const auto name=object.parameters.at("preset").get<std::string>();
        if(name!="singleFluidLaminar" && name!="singleFluidSST")
            throw std::runtime_error(object.source+": unknown composition preset "+name);
        const bool sst=name=="singleFluidSST";
        const auto source=object.source+" ["+name+"]";
        merge(merge,numerical,P{{"time",{{"default",sst?"classicalRK4":"forwardEuler"}}},
            {"terms",{{"convection","weno5Rusanov"},{"diffusion","central2Explicit"}}}},"solver/numerics",source,false);
        const auto add=[&](const std::string& id,const std::string& type,P parameters) {
            auto& value=defaults[role(type)];
            if(value.type.empty())value={id,"solver",type};
            if(type=="turbulence")value.category="models";
            merge(merge,value.parameters,parameters,"object/"+role(type),source,false);
            value.source=source;
        };
        add("equations","equationSystem",{{"preset","compressible"}});
        add("state","stateRegistry",{{"use",P::array({"rho","rhoU","rhoE"})}});
        add("equationRegistry","equationRegistry",{{"use",P::array({"Continuity","Momentum","Energy"})}});
        add("execution","algorithmRegistry",{{"Explicit",P::object()}});
        add("turbulence","turbulence",sst?P{{"simulationType","RAS"},{"RAS",{{"model","kOmegaSST"}}},{"transportMethod","stageConservative"}}:P{{"simulationType","laminar"}});
    }
    const auto numericsSource=input.selectionOrigins.find("solver/numerics");
    merge(merge,numerical,input.numerics,"solver/numerics",numericsSource==input.selectionOrigins.end()?"registered numerics":numericsSource->second,true);
    result.numerics=std::move(numerical);
    result.objects.clear();
    std::map<std::string,std::string> explicitRoles;
    for(auto object:explicitObjects) {
        const auto key=role(object.type);
        if(defaults.count(key)) {
            const auto [prior,inserted]=explicitRoles.emplace(key,object.source);
            if(!inserted)throw std::runtime_error("Conflicting explicit target "+key+": "+prior->second+" and "+object.source);
            auto parameters=defaults.at(key).parameters;
            // Explicit algorithm/model choices replace alternative branches, not merge two algorithms.
            if(key=="algorithmregistry" || (key=="turbulence" && object.parameters.contains("simulationType"))) {
                for(const auto& item:parameters.items())if(!object.parameters.contains(item.key())) {
                    const auto prefix="object/"+key+"/"+item.key();
                    for(auto it=origins.begin();it!=origins.end();) {
                        if(it->first==prefix || it->first.rfind(prefix+"/",0)==0)it=origins.erase(it);else ++it;
                    }
                }
                parameters=P::object();
            }
            merge(merge,parameters,object.parameters,"object/"+key,object.source,true);
            object.parameters=std::move(parameters);
        } else {
            P empty;merge(merge,empty,object.parameters,"object/"+key,object.source,true);
        }
        result.objects.push_back(std::move(object));
    }
    for(auto& [key,object]:defaults)if(!explicitRoles.count(key))result.objects.push_back(std::move(object));
    origins["solver/numerics"]=numericsSource==input.selectionOrigins.end()?"registered numerics":numericsSource->second;
    if(!result.numerics.contains("time"))origins["solver/numerics/time/default"]="builtin NumericsConfig default; no explicit time selection";
    result.selectionOrigins=std::move(origins);
    return result;
}
} // namespace
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
    auto builtins=expandPresets(model);
    const auto instances = equationInstances(builtins);
    const auto objects=builtins.objects;
    std::vector<Model::ObjectDescriptor> extensions;
    builtins.objects.clear();
    for(const auto& object:objects) {
        if(objectRegistry().contains(object.type))extensions.push_back(object);
        else builtins.objects.push_back(object);
    }
    auto result=CaseAdapter(model.caseDir).build(builtins);
    for(const auto& object:extensions)objectRegistry().create(object,result);
    auto authored=model;authored.selectionOrigins=builtins.selectionOrigins;
    result.modelDescription=std::make_shared<const Model::Description>(std::move(authored));
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
