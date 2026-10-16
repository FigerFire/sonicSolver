#include "SF_IORegistry.h"
#include <algorithm>
#include <cctype>
namespace SF::IO {
namespace {
std::string normalized(std::string value) {
    std::transform(value.begin(),value.end(),value.begin(),[](unsigned char c){return static_cast<char>(std::tolower(c));});
    return value;
}
class FunctionHandler final : public IOHandler {
public:
    explicit FunctionHandler(IORegistry::Loader loader):loader_(std::move(loader)) {}
    void load(LoadContext& context,const Serialization::SonicDocument& document) const override {
        loader_(context,document);
    }
private:
    IORegistry::Loader loader_;
};
}
void IORegistry::add(const std::string& object,const std::string& type,Loader loader) {
    if(!loader)throw std::runtime_error("Cannot register an empty IO handler");
    const Key key{normalized(object),normalized(type)};
    if(!handlers_.emplace(key,std::make_shared<FunctionHandler>(std::move(loader))).second)
        throw std::runtime_error("Duplicate IO handler: "+object+"/"+type);
}
void IORegistry::addDefaultFile(const std::string& object,const std::string& entryName,const std::string& file) {
    const Key key{normalized(object),normalized(entryName)};
    if(!defaultFiles_.emplace(key,file).second)throw std::runtime_error("Duplicate default SonicFile: "+object+"/"+entryName);
}
void IORegistry::reserveEntry(const std::string& object,const std::string& name,const std::string& role,bool singleton) {
    if(!entryRules_.emplace(Key{normalized(object),normalized(name)},EntryRule{name,role,singleton}).second)
        throw std::runtime_error("Duplicate reserved registry entry: "+object+"/"+name);
}
const IOHandler& IORegistry::find(const std::string& object,const std::string& type) const {
    const auto name=normalized(object),kind=normalized(type);
    auto it=handlers_.find({name,kind});
    if(it==handlers_.end())it=handlers_.find({name,"*"});
    if(it==handlers_.end())throw std::runtime_error("No IO handler registered for object "+object+" type "+type);
    return *it->second;
}
RegistryEntry IORegistry::parseEntry(const std::string& object,const std::string& name,const Model::Parameters& value,const std::string& source) const {
    RegistryEntry result;result.name=name;
    if(value.is_string()) {
        // A scalar is a compact built-in selection: multiPhase: eulerianEulerian.
        // This is still an explicit registration, never a search for unregistered files.
        result.type=name;result.selection=value.get<std::string>();
    } else {
        Model::Schema{{{"type","string",true},{"file","string"},{"selection","string"}},false}.validate(value,source+"/"+name);
        result.type=value.at("type").get<std::string>();result.file=value.value("file",std::string());result.selection=value.value("selection",std::string());
    }
    if(name.empty() || result.type.empty() || normalized(name)=="registry" || normalized(result.type)=="registry")throw std::runtime_error(source+"/"+name+": invalid/reserved registry identity");
    const auto rule=entryRules_.find({normalized(object),normalized(name)});
    if(rule!=entryRules_.end() && rule->second.name!=name)
        throw std::runtime_error(source+"/"+name+": reserved registry key requires exact identity "+rule->second.name);
    if(rule!=entryRules_.end() && rule->second.role!=result.type)
        throw std::runtime_error(source+"/"+name+": reserved registry role requires type "+rule->second.role+", got "+result.type);
    for(const auto& item:entryRules_)if(item.second.singleton && normalized(item.second.role)==normalized(result.type)) {
        if(item.second.role!=result.type)throw std::runtime_error(source+"/"+name+": reserved role requires exact type "+item.second.role);
        if(item.first.object!=normalized(object))throw std::runtime_error(source+"/"+name+": reserved role "+result.type+" belongs to object "+item.first.object);
        if(item.first.type!=normalized(name))throw std::runtime_error(source+"/"+name+": reserved role "+result.type+" requires registry key "+item.first.type);
    }
    return result;
}
std::string IORegistry::resolveFile(const std::string& object,const RegistryEntry& entry,const std::string& source) const {
    if(!entry.file.empty())return entry.file;
    const auto it=defaultFiles_.find({normalized(object),normalized(entry.name)});
    if(it!=defaultFiles_.end())return it->second;
    throw std::runtime_error(source+"/"+entry.name+": custom registry entry requires file");
}
} // namespace SF::IO
