#include "core/model/SF_model.h"
namespace SF::Model {
void Schema::validate(const Parameters& value,const std::string& context) const {
    if(!value.is_object()) throw std::runtime_error(context+": parameters must be a mapping");
    for(const auto& rule:parameters) {
        auto it=value.find(rule.name);
        if(it==value.end()) { if(rule.required) throw std::runtime_error(context+": missing "+rule.name); continue; }
        const bool valid=rule.type=="any" || (rule.type=="number"&&it->is_number())
            || (rule.type=="integer"&&it->is_number_integer()) || (rule.type=="string"&&it->is_string())
            || (rule.type=="boolean"&&it->is_boolean()) || (rule.type=="array"&&it->is_array())
            || (rule.type=="object"&&it->is_object());
        if(!valid) throw std::runtime_error(context+"/"+rule.name+": expected "+rule.type);
    }
    if(!allowAdditional) for(auto it=value.begin();it!=value.end();++it) {
        bool known=false; for(const auto& rule:parameters) known|=rule.name==it.key();
        if(!known) throw std::runtime_error(context+": unknown parameter "+it.key());
    }
}
void FieldRegistry::add(FieldView field) {
    const auto name=field.descriptor.name;
    if(name.empty()||!field.read) throw std::runtime_error("Field registration requires name and reader");
    if(!fields_.emplace(name,std::move(field)).second) throw std::runtime_error("Duplicate field: "+name);
}
FieldView& FieldRegistry::lookup(const std::string& name) {
    auto it=fields_.find(name);if(it==fields_.end())throw std::runtime_error("Unknown field: "+name);
    return it->second;
}
}
