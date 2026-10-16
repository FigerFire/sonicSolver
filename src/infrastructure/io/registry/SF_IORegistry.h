#pragma once
/// @file SF_IORegistry.h
/// @brief SonicFile 的 object/注册角色分派表。它只管理文件语义，不包含物理模型实现。
#include "core/model/SF_model.h"
#include "infrastructure/io/serialization/SF_serialization.h"
#include <filesystem>
#include <functional>
#include <memory>
#include <optional>
#include <string>

namespace SF::IO {
class IORegistry;
/// @brief name 是文件身份，type 是内部角色；二者独立校验，file 只由显式注册解析。
struct RegistryEntry {
    std::string name, type, file, selection;
};
/// @brief 一次文件加载的上下文。loadFile 负责再次读取标题并分派给对应 handler。
struct LoadContext {
    Model::Description& model;
    const IORegistry& registry;
    std::filesystem::path caseDirectory;
    std::string entryName, entryType, selection;
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
/// @brief object/内部角色到 handler 的注册表；* 是该 object 的通用模块 handler。
class IORegistry {
public:
    using Loader=std::function<void(LoadContext&,const Serialization::SonicDocument&)>;
    void add(const std::string& object,const std::string& type,Loader loader);
    /// @brief 给内置 registry 条目登记默认模块位置；自定义条目必须显式声明 file。
    void addDefaultFile(const std::string& object,const std::string& entryName,const std::string& file);
    /// @brief 保留注册键及其内部角色；singleton 角色不能换名绕过唯一 authority。
    void reserveEntry(const std::string& object,const std::string& name,const std::string& role,bool singleton=false);
    const IOHandler& find(const std::string& object,const std::string& type) const;
    RegistryEntry parseEntry(const std::string& object,const std::string& name,const Model::Parameters& value,
                             const std::string& source) const;
    std::string resolveFile(const std::string& object,const RegistryEntry& entry,
                            const std::string& source) const;
private:
    struct Key { std::string object,type; bool operator<(const Key& other) const {
        return object<other.object||(object==other.object&&type<other.type);
    }};
    std::map<Key,std::shared_ptr<const IOHandler>> handlers_;
    std::map<Key,std::string> defaultFiles_;
    struct EntryRule {std::string name,role;bool singleton;};
    std::map<Key,EntryRule> entryRules_;
};
} // namespace SF::IO
