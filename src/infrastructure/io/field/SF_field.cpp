#include "SF_field.h"
#include "infrastructure/io/serialization/SF_serialization.h"
#include <cmath>
#include <set>
namespace SF::FieldIO {
namespace {
Model::ValueType valueType(const std::string& type,const std::string& name) {
    if(type=="scalar")return Model::ValueType::Scalar;
    if(type=="vector")return Model::ValueType::Vector;
    if(type=="symmetricTensor")return Model::ValueType::SymmetricTensor;
    if(type=="tensor")return Model::ValueType::Tensor;
    throw std::runtime_error("field/"+name+": unknown value type "+type);
}
void validateDescriptor(Model::FieldDescriptor& f,const std::string& name) {
    if(f.storage!="primary"&&f.storage!="solveVariable"&&f.storage!="derived")throw std::runtime_error("field/"+name+": invalid storage role");
    if(!f.dimensions.is_array()||f.dimensions.size()!=7)throw std::runtime_error("field/"+name+": dimensions must contain 7 values");
    for(const auto& x:f.dimensions)if(!x.is_number())throw std::runtime_error("field/"+name+": dimensions must be numeric");
}
}
Model::FieldDescriptor decodeRegistry(const std::string& name,const Model::Parameters& v) {
    if(!v.is_object())throw std::runtime_error("field/"+name+": registry entry must be a mapping");
    Model::Schema schema{{{"type","string",true},{"domain","string"},{"location","string"},{"dimensions","array"},
        {"storage","string"},{"output","boolean"},{"outputName","string"},{"restart","boolean"},
        {"new","boolean"}},false};
    schema.validate(v,"field/"+name);
    Model::FieldDescriptor f;f.name=name;f.symbol=name;f.type=valueType(v.at("type").get<std::string>(),name);
    f.domain=v.value("domain",f.domain);f.location=v.value("location",f.location);
    f.dimensions=v.value("dimensions",f.dimensions);f.storage=v.value("storage",f.storage);
    f.output=v.value("output",true);f.restart=v.value("restart",false);f.outputName=v.value("outputName",name);
    validateDescriptor(f,name);return f;
}
Model::FieldDescriptor decode(const std::string& name,const Model::Parameters& v) {
    Model::FieldDescriptor f;f.name=name;f.symbol=name;Model::Parameters initial=v;
    if(v.is_object()&&v.contains("type")) {
        auto descriptor=v;
        descriptor.erase("default");descriptor.erase("value");descriptor.erase("initial");descriptor.erase("sets");
        f=decodeRegistry(name,descriptor);
        if(v.contains("default"))initial=v.at("default");else if(v.contains("value"))initial=v.at("value");else if(v.contains("initial"))initial=v.at("initial");else initial=nullptr;
        const std::set<std::string> reserved{"type","domain","location","dimensions","storage","output","outputName","restart","default","value","initial","sets"};
        f.sets=v.value("sets",Model::Parameters::object());
        for(auto it=v.begin();it!=v.end();++it)if(!reserved.count(it.key()))f.sets[it.key()]=it.value().is_object()&&it.value().contains("type")?it.value():Model::Parameters{{"type","fixedValue"},{"value",it.value()}};
    } else {
        f.outputName=name;f.output=true;f.restart=false;
        f.type=initial.is_array()?(initial.size()==3?Model::ValueType::Vector:initial.size()==6?Model::ValueType::SymmetricTensor:Model::ValueType::Tensor):Model::ValueType::Scalar;
    }
    if(initial.is_null())throw std::runtime_error("field/"+name+": missing initial value");
    f.initial=initial;
    const int n=Model::components(f.type);
    if(n==1&&!f.initial.is_number())throw std::runtime_error("field/"+name+": initial must be numeric");
    if(n!=1&&(!f.initial.is_array()||f.initial.size()!=std::size_t(n)))throw std::runtime_error("field/"+name+": initial component count mismatch");
    if(n!=1)for(const auto& x:f.initial)if(!x.is_number())throw std::runtime_error("field/"+name+": nonnumeric component");
    return f;
}
Model::Parameters encode(const Model::FieldDescriptor& f) {
    const char* types[]={"scalar","vector","symmetricTensor","tensor"};
    return {{"type",types[static_cast<int>(f.type)]},{"default",f.initial},{"output",f.output},{"restart",f.restart}};
}
Model::Parameters encodeRegistry(const Model::FieldDescriptor& f) {
    const char* types[]={"scalar","vector","symmetricTensor","tensor"};
    Model::Parameters value={{"type",types[static_cast<int>(f.type)]},{"domain",f.domain},
        {"location",f.location},{"dimensions",f.dimensions},{"storage",f.storage}};
    if(!f.output) value["output"]=false;
    if(f.outputName!=f.name&&!f.outputName.empty()) value["outputName"]=f.outputName;
    if(f.restart) value["restart"]=true;
    return value;
}
void write(const Model::FieldRegistry& registry,const std::string& file,bool restartOnly) {
    auto data=Model::Parameters::object();
    for(const auto& [name,f]:registry.fields()) {
        if(restartOnly&&!f.descriptor.restart)continue;
        auto values=Model::Parameters::array();
        for(std::size_t i=0;i<f.size;++i)for(int c=0;c<Model::components(f.descriptor.type);++c) {
            double v=f.read(i,c);if(!std::isfinite(v))throw std::runtime_error("Nonfinite field checkpoint: "+name);
            values.push_back(v);
        }
        data[name]={{"descriptor",encode(f.descriptor)},{"size",f.size},{"values",values}};
    }
    Serialization::write({{"formatVersion",1},{"fields",data}},file);
}
void read(Model::FieldRegistry& registry,const std::string& file,bool restartOnly) {
    const auto document=Serialization::read(file);
    if(document.at("formatVersion")!=1)throw std::runtime_error(file+": unsupported field version");
    const auto& data=document.at("fields");
    // Validate the whole checkpoint before modifying any registered state.
    for(const auto& [name,f]:registry.fields()) {
        if(restartOnly&&!f.descriptor.restart)continue;
        if(!data.contains(name)||!f.write)throw std::runtime_error(file+": missing/unwritable field "+name);
        const auto& item=data.at(name);const auto d=decode(name,item.at("descriptor"));
        if(d.type!=f.descriptor.type
            ||item.at("size")!=f.size||item.at("values").size()!=f.size*Model::components(d.type))
            throw std::runtime_error(file+": incompatible field "+name);
        for(const auto& v:item.at("values"))if(!v.is_number()||!std::isfinite(v.get<double>()))throw std::runtime_error(file+": invalid array "+name);
    }
    for(const auto& [name,f]:registry.fields()) {
        if(restartOnly&&!f.descriptor.restart)continue;const auto& values=data.at(name).at("values");
        const int n=Model::components(f.descriptor.type);
        for(std::size_t i=0;i<f.size;++i)for(int c=0;c<n;++c)f.write(i,c,values.at(i*n+c).get<double>());
    }
}
}
