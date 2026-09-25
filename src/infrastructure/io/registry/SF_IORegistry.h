#pragma once
/// @file SF_IORegistry.h
/// @brief SonicFile 的 object/type 分派表。它只管理文件语义，不包含物理模型实现。
#include "core/model/SF_model.h"
#include "infrastructure/io/serialization/SF_serialization.h"
#include <filesystem>
#include <functional>
#include <memory>
#include <optional>
#include <string>

namespace SF::IO {
class IORegistry;
/// @brief Registry 条目的通用引用；type 选择处理器，file 指向可替换的模块配置。
struct RegistryEntry {
    std::string name, type, file, selection;
};
/// @brief 一次文件加载的上下文。loadFile 负责再次读取标题并分派给对应 handler。
struct LoadContext {
    Model::Description& model;
    const IORegistry& registry;
    std::filesystem::path caseDirectory;
    std::string entryName, selection;
    std::function<void(const std::string& relative,const std::string& entryName,
                       const std::string& expectedObject,const std::string& expectedType,
                       const std::string& selection)> loadFile;
};
/// @brief 任意 SonicFile handler 的小接口。
class IOHandler {
public:
    virtual ~IOHandler()=default;
    virtual void load(LoadContext& context,const Serialization::SonicDocument& document) const=0;
};
/// @brief object/type 到 handler 的注册表；* 是该 object 的通用模块 handler。
class IORegistry {
public:
    using Loader=std::function<void(LoadContext&,const Serialization::SonicDocument&)>;
    void add(const std::string& object,const std::string& type,Loader loader);
    /// @brief 给内置 registry 条目登记默认模块位置；自定义条目必须显式声明 file。
    void addDefaultFile(const std::string& object,const std::string& entryName,const std::string& file);
    const IOHandler& find(const std::string& object,const std::string& type) const;
    RegistryEntry parseEntry(const std::string& name,const Model::Parameters& value,
                             const std::string& source) const;
    std::string resolveFile(const std::string& object,const RegistryEntry& entry,
                            const std::string& source) const;
private:
    struct Key { std::string object,type; bool operator<(const Key& other) const {
        return object<other.object||(object==other.object&&type<other.type);
    }};
    std::map<Key,std::shared_ptr<const IOHandler>> handlers_;
    std::map<Key,std::string> defaultFiles_;
};
} // namespace SF::IO
