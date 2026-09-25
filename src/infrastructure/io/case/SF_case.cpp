#include "SF_case.h"
#include "infrastructure/io/field/SF_field.h"
#include "infrastructure/io/registry/SF_IORegistry.h"
#include "infrastructure/io/serialization/SF_serialization.h"
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <optional>
#include <set>

namespace SF::CaseIO {
namespace {
using P=Model::Parameters;
using Sonic=Serialization::SonicDocument;

std::string lower(std::string value) {
    for(char& c:value)c=static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return value;
}
bool same(const std::string& a,const std::string& b) { return b.empty()||lower(a)==lower(b); }
P compactBoundary(const P& value) {
    if(value.is_string())return {{"type",value.get<std::string>()}};
    if(!value.is_object())throw std::runtime_error("Boundary value must be a condition name or mapping");
    if(value.contains("type"))return value;
    if(value.size()!=1)throw std::runtime_error("Boundary condition must contain exactly one condition type");
    const auto item=value.begin();
    return item.value().is_null()?P{{"type",item.key()}}:P{{"type",item.key()},{"value",item.value()}};
}
P compactBoundaries(const P& value) {
    if(!value.is_object())throw std::runtime_error("Field boundaries must be a mapping");
    P result=P::object();
    for(auto it=value.begin();it!=value.end();++it)result[it.key()]=compactBoundary(it.value());
    return result;
}
P simpleBoundary(const P& value) {
    if(!value.is_object()||!value.contains("type"))return value;
    if(!value.contains("value"))return value.at("type");
    return {{value.at("type").get<std::string>(),value.at("value")}};
}
P simpleBoundaries(const P& value) {
    if(!value.is_object())return value;
    P result=P::object();
    for(auto field=value.begin();field!=value.end();++field) {
        if(field.key()=="default") { result[field.key()]=simpleBoundary(field.value());continue; }
        P patches=P::object();
        for(auto patch=field.value().begin();patch!=field.value().end();++patch)
            patches[patch.key()]=simpleBoundary(patch.value());
        result[field.key()]=std::move(patches);
    }
    return result;
}
void addObject(Model::Description& model,Model::ObjectDescriptor object,const std::string& source) {
    const auto duplicate=std::find_if(model.objects.begin(),model.objects.end(),[&](const auto& current) {
        return current.name==object.name&&current.category==object.category;
    });
    if(duplicate!=model.objects.end())throw std::runtime_error(source+": duplicate "+object.category+" object "+object.name);
    model.objects.push_back(std::move(object));
}
std::string safeFileName(std::string value) {
    for(char& c:value)if(!std::isalnum(static_cast<unsigned char>(c))&&c!='-'&&c!='_')c='_';
    return value.empty()?"module":value;
}
std::string semanticType(const std::string& type) {
    static const std::map<std::string,std::string> legacy{
        {"BuiltinEquationSystem","equationSystem"},{"ImmersedBoundary","IBM"},
        {"PhaseSystem","multiPhase"},{"Thermophysical","thermophysical"},
        {"PhaseChange","phaseChange"},{"Turbulence","turbulence"},
        {"Gravity","gravity"},{"WallHeat","wallHeat"},{"File","geometry"}};
    const auto found=legacy.find(type);return found==legacy.end()?type:found->second;
}
bool solverModule(const Model::ObjectDescriptor& object) {
    const auto type=semanticType(object.type);
    return object.category=="solver"||type=="IBM"||type=="ILW"||type=="equationSystem";
}

std::optional<Model::FieldDescriptor> builtinField(
        const std::string& name) {
    using Type = Model::ValueType;
    struct Entry { const char* name; Type type; };
    static const Entry catalog[] = {
        {"rho",Type::Scalar},{"U",Type::Vector},{"p",Type::Scalar},
        {"T",Type::Scalar},{"h",Type::Scalar},{"e",Type::Scalar},
        {"alpha",Type::Scalar},{"k",Type::Scalar},{"omega",Type::Scalar},
        {"epsilon",Type::Scalar}
    };
    for (const auto& entry : catalog) {
        if (name == entry.name) {
            Model::FieldDescriptor result;
            result.name = name;
            result.symbol = name;
            result.outputName = name;
            result.type = entry.type;
            return result;
        }
    }
    return std::nullopt;
}

IO::IORegistry nativeRegistry() {
    IO::IORegistry registry;
    registry.add("fields","registry",[](IO::LoadContext& context,const Sonic& document) {
        if(!document.body.is_object())throw std::runtime_error(document.source+": fields registry body must be a mapping");
        const auto defaultValue=document.body.contains("default")?document.body.at("default"):P();
        std::set<std::string> symbols;
        for(auto it=document.body.begin();it!=document.body.end();++it) {
            if(it.key()=="default")continue;
            const auto& value=it.value().is_null()&&!defaultValue.is_null()?defaultValue:it.value();
            Model::FieldDescriptor field;
            if (value.is_object() && !value.contains("type")) {
                const auto builtin = builtinField(it.key());
                if (!builtin) {
                    throw std::runtime_error(
                        document.source+": custom field '"+it.key()
                        +"' must declare new: true and its type.");
                }
                static const std::set<std::string> allowed{
                    "output","outputName","restart"};
                for (auto property = value.begin();property != value.end();++property) {
                    if (!allowed.count(property.key())) {
                        throw std::runtime_error(
                            document.source+": built-in field '"+it.key()
                            +"' only accepts IO metadata overrides.");
                    }
                }
                field = *builtin;
                field.output = value.value("output",field.output);
                field.outputName = value.value("outputName",field.outputName);
                field.restart = value.value("restart",field.restart);
            } else {
                if (value.is_object() && value.contains("new")
                    && !value.at("new").get<bool>()) {
                    throw std::runtime_error(
                        document.source+": custom field '"+it.key()
                        +"' must set new: true.");
                }
                field=FieldIO::decodeRegistry(it.key(),value);
            }
            if(!symbols.insert(field.symbol).second)throw std::runtime_error(document.source+": duplicate field symbol "+field.symbol);
            context.model.fields.push_back(std::move(field));
        }
    });
    registry.add("fields","internalField",[](IO::LoadContext& context,const Sonic& document) {
        if(!document.body.is_object())throw std::runtime_error(document.source+": internalField body must be a mapping");
        for(auto it=document.body.begin();it!=document.body.end();++it) {
            auto found=std::find_if(context.model.fields.begin(),context.model.fields.end(),[&](const auto& field){return field.name==it.key();});
            if(found==context.model.fields.end())throw std::runtime_error(document.source+": unknown field "+it.key());
            const auto& value=it.value();
            if(value.is_object()&&value.contains("uniform"))found->initial=value.at("uniform");
            else if(value.is_object()&&value.contains("default")) {
                found->initial=value.at("default");found->sets=P::object();
                for(auto set=value.begin();set!=value.end();++set)if(set.key()!="default")
                    found->sets[set.key()]=set.value().is_object()&&set.value().contains("type")?set.value():P{{"type","fixedValue"},{"value",set.value()}};
            } else found->initial=value;
        }
    });
    registry.add("fields","boundary",[](IO::LoadContext& context,const Sonic& document) {
        if(!document.body.is_object())throw std::runtime_error(document.source+": boundary body must be a mapping");
        const auto defaultValue=document.body.contains("default")?compactBoundary(document.body.at("default")):P();
        for(auto it=document.body.begin();it!=document.body.end();++it) {
            if(it.key()=="default")continue;
            auto found=std::find_if(context.model.fields.begin(),context.model.fields.end(),[&](const auto& field){return field.name==it.key();});
            if(found==context.model.fields.end())throw std::runtime_error(document.source+": boundary references unknown field "+it.key());
            if(!found->boundaries.empty())throw std::runtime_error(document.source+": duplicate boundaries for "+it.key());
            found->boundaries=it.value().is_null()?P::object():compactBoundaries(it.value());
            if(!defaultValue.is_null()&&!found->boundaries.contains("default"))found->boundaries["default"]=defaultValue;
        }
    });
    const auto loadRegistry=[](IO::LoadContext& context,const Sonic& document) {
        if(!document.body.is_object())throw std::runtime_error(document.source+": registry body must be a mapping");
        for(auto it=document.body.begin();it!=document.body.end();++it) {
            const auto entry=context.registry.parseEntry(it.key(),it.value(),document.source);
            const auto file=context.registry.resolveFile(document.object,entry,document.source);
            context.loadFile(file,entry.name,document.object,entry.type,entry.selection);
        }
    };
    registry.add("solver","registry",loadRegistry);
    registry.add("models","registry",loadRegistry);
    registry.add("solver","runtime",[](IO::LoadContext& context,const Sonic& document) {
        if(!context.model.runtime.empty())throw std::runtime_error(document.source+": runtime is already loaded");
        context.model.runtime=document.body;
    });
    registry.add("solver","numerics",[](IO::LoadContext& context,const Sonic& document) {
        if(!context.model.numerics.empty())throw std::runtime_error(document.source+": numerics is already loaded");
        context.model.numerics=document.body;
    });
    registry.add("solver","algorithm",[](IO::LoadContext& context,const Sonic& document) {
        if(!context.model.solver.empty())throw std::runtime_error(document.source+": algorithm is already loaded");
        context.model.solver=document.body;
    });
    registry.add("equations","registry",[](IO::LoadContext& context,const Sonic& document) {
        addObject(context.model,{"equations","equations","equationRegistry",{},document.body,{}},document.source);
    });
    registry.add("algorithms","registry",[](IO::LoadContext& context,const Sonic& document) {
        addObject(context.model,{"algorithms","algorithms","algorithmRegistry",{},document.body,{}},document.source);
    });
    registry.add("solver","*",[](IO::LoadContext& context,const Sonic& document) {
        if(context.entryName.empty())throw std::runtime_error(document.source+": a solver module requires a registry name");
        auto parameters=document.body;const auto expression=parameters.value("expression",std::string());
        if(!expression.empty()) {
            parameters.erase("expression");
            if(parameters.contains("parameters"))parameters=parameters.at("parameters");
        }
        addObject(context.model,{context.entryName,"solver",document.type,expression,parameters,context.selection},document.source);
    });
    registry.add("models","*",[](IO::LoadContext& context,const Sonic& document) {
        if(context.entryName.empty())throw std::runtime_error(document.source+": a model module requires a registry name");
        auto parameters=document.body;const auto expression=parameters.value("expression",std::string());
        if(!expression.empty()) {
            parameters.erase("expression");
            if(parameters.contains("parameters"))parameters=parameters.at("parameters");
        }
        addObject(context.model,{context.entryName,"models",document.type,expression,parameters,context.selection},document.source);
    });
    registry.add("mesh","registry",[](IO::LoadContext& context,const Sonic& document) {
        if(!context.model.mesh.empty())throw std::runtime_error(document.source+": mesh is already loaded");
        context.model.mesh=document.body;
    });
    // Declarations only: the generic registry resolver needs no IBM/ILW branch.
    registry.addDefaultFile("solver","IBM","models/IBM.yaml");
    registry.addDefaultFile("solver","ILW","models/ILW.yaml");
    registry.addDefaultFile("models","multiPhase","models/multiPhase.yaml");
    registry.addDefaultFile("models","turbulence","models/turbulence.yaml");
    registry.addDefaultFile("models","thermoDynamics","models/thermoDynamics.yaml");
    return registry;
}
} // namespace

Model::Description read(const std::string& path) {
    namespace fs=std::filesystem;
    const fs::path directory=fs::is_directory(path)?fs::path(path):fs::path(path).parent_path();
    Model::Description model;model.caseDir=directory.string();model.name=directory.filename().string();
    const auto registry=nativeRegistry();
    std::function<void(const std::string&,const std::string&,const std::string&,const std::string&,const std::string&)> dispatch;
    dispatch=[&](const std::string& relative,const std::string& entryName,const std::string& expectedObject,const std::string& expectedType,const std::string& selection) {
        const auto full=(directory/relative).lexically_normal();
        const auto document=Serialization::readSonicFile(full.string());
        const auto shown=fs::relative(full,directory).string();
        model.comments[shown]=document.comments;
        if(!same(document.object,expectedObject)||!same(document.type,expectedType))
            throw std::runtime_error(shown+": expected object "+expectedObject+" type "+expectedType+
                                     ", got object "+document.object+" type "+document.type);
        IO::LoadContext context{model,registry,directory,entryName,selection,dispatch};
        registry.find(document.object,document.type).load(context,document);
    };
    const auto meta=Serialization::readSonicFile((directory/"case.yaml").string());
    model.comments["case.yaml"]=meta.comments;
    if(!same(meta.object,"case")||!same(meta.type,"registry"))throw std::runtime_error("case.yaml: expected object case type registry");
    Model::Schema{{{"formatVersion","integer",true},{"name","string",true}},false}.validate(meta.body,"case.yaml");
    model.formatVersion=meta.body.at("formatVersion");model.name=meta.body.at("name");
    if(model.formatVersion!=1)throw std::runtime_error("case.yaml: unsupported formatVersion "+std::to_string(model.formatVersion));
    dispatch("mesh/mesh.yaml","","mesh","registry","");
    dispatch("fields/fields.yaml","","fields","registry","");
    dispatch("fields/internalField.yaml","","fields","internalField","");
    dispatch("fields/boundaries.yaml","","fields","boundary","");
    dispatch("solvers/solvers.yaml","","solver","registry","");
    if (fs::exists(directory/"equations/equations.yaml")) {
        dispatch("equations/equations.yaml","","equations","registry","");
    }
    if (fs::exists(directory/"algorithms/algorithms.yaml")) {
        dispatch("algorithms/algorithms.yaml","","algorithms","registry","");
    }
    dispatch("models/models.yaml","","models","registry","");
    return model;
}

void write(const Model::Description& model,const std::string& directory) {
    const auto dir=std::filesystem::path(directory);
    auto writeFile=[&](const std::string& relative,const std::string& object,const std::string& type,
                       const P& body,bool compactConditions=false) {
        const auto found=model.comments.find(relative);
        Serialization::writeSonicFile({object,type,(dir/relative).string(),body,
            found==model.comments.end()?std::vector<std::string>{}:found->second},(dir/relative).string(),compactConditions);
    };
    writeFile("case.yaml","case","registry",{{"formatVersion",model.formatVersion},{"name",model.name}});
    writeFile("mesh/mesh.yaml","mesh","registry",model.mesh);
    P fieldRegistry=P::object();
    for(const auto& field:model.fields)fieldRegistry[field.name]=FieldIO::encodeRegistry(field);
    writeFile("fields/fields.yaml","fields","registry",fieldRegistry);
    P internal=P::object();
    for(const auto& field:model.fields) {
        if(field.initial.is_null())continue;
        if(field.sets.empty())internal[field.name]={{"uniform",field.initial}};
        else {
            P value={{"default",field.initial}};
            for(auto set=field.sets.begin();set!=field.sets.end();++set)value[set.key()]=set.value().value("value",set.value());
            internal[field.name]=std::move(value);
        }
    }
    writeFile("fields/internalField.yaml","fields","internalField",internal);
    P boundaries=P::object();for(const auto& field:model.fields)boundaries[field.name]=simpleBoundaries(field.boundaries);
    writeFile("fields/boundaries.yaml","fields","boundary",boundaries,true);

    P solverRegistry=P::object();
    solverRegistry["runtime"]={{"type","runtime"},{"file","solvers/runtime.yaml"}};
    solverRegistry["numerics"]={{"type","numerics"},{"file","solvers/numerics.yaml"}};
    solverRegistry["algorithm"]={{"type","algorithm"},{"file","solvers/algorithm.yaml"}};
    writeFile("solvers/runtime.yaml","solver","runtime",model.runtime);
    writeFile("solvers/numerics.yaml","solver","numerics",model.numerics);
    writeFile("solvers/algorithm.yaml","solver","algorithm",model.solver);
    P modelRegistry=P::object();
    std::set<std::string> usedSolverNames,usedModelNames;
    for(const auto& object:model.objects) {
        const auto type=semanticType(object.type);
        const bool solver=solverModule(object);
        auto& used=solver?usedSolverNames:usedModelNames;
        if(!used.insert(object.name).second)throw std::runtime_error("Duplicate registry object: "+object.name);
        const std::string relative="models/"+safeFileName(object.name)+".yaml";
        P body=object.parameters.is_null()?P::object():object.parameters;
        if(!object.expression.empty())body={{"expression",object.expression},{"parameters",body}};
        writeFile(relative,solver?"solver":"models",type,body);
        P entry={{"type",type},{"file",relative}};
        if(solver)solverRegistry[object.name]=std::move(entry);else modelRegistry[object.name]=std::move(entry);
    }
    writeFile("solvers/solvers.yaml","solver","registry",solverRegistry);
    writeFile("models/models.yaml","models","registry",modelRegistry);
}
} // namespace SF::CaseIO
