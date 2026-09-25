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
const IOHandler& IORegistry::find(const std::string& object,const std::string& type) const {
    const auto name=normalized(object),kind=normalized(type);
    auto it=handlers_.find({name,kind});
    if(it==handlers_.end())it=handlers_.find({name,"*"});
    if(it==handlers_.end())throw std::runtime_error("No IO handler registered for object "+object+" type "+type);
    return *it->second;
}
RegistryEntry IORegistry::parseEntry(const std::string& name,const Model::Parameters& value,const std::string& source) const {
    RegistryEntry result;result.name=name;
    if(value.is_string()) {
        // A scalar is a compact built-in selection: multiPhase: eulerianEulerian.
        // The registry key remains the module type and therefore still matches the file title.
        result.type=name;result.selection=value.get<std::string>();
        return result;
    }
    Model::Schema{{{"type","string",true},{"file","string"}},false}.validate(value,source+"/"+name);
    result.type=value.at("type").get<std::string>();result.file=value.value("file",std::string());
    return result;
}
std::string IORegistry::resolveFile(const std::string& object,const RegistryEntry& entry,const std::string& source) const {
    if(!entry.file.empty())return entry.file;
    const auto it=defaultFiles_.find({normalized(object),normalized(entry.name)});
    if(it!=defaultFiles_.end())return it->second;
    throw std::runtime_error(source+"/"+entry.name+": custom registry entry requires file");
}
} // namespace SF::IO
