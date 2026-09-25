#pragma once
/// @file SF_model.h
/// @brief 前端共享的通用模型描述；不依赖具体物理、求解器或磁盘格式。
#include <nlohmann/json.hpp>
#include <map>
#include <string>
#include <vector>
#include <functional>
#include <stdexcept>
#include <optional>
namespace SF::Model {
using Parameters = nlohmann::ordered_json;
/// @brief 字段基础值类型，与内存布局独立。
enum class ValueType { Scalar, Vector, SymmetricTensor, Tensor };
inline int components(ValueType t) {
    switch(t) { case ValueType::Scalar:return 1; case ValueType::Vector:return 3;
        case ValueType::SymmetricTensor:return 6; case ValueType::Tensor:return 9; }
    throw std::runtime_error("Invalid field value type");
}
/// @brief 通用字段元数据。registry 文件声明语义，internalField 文件提供数值。
struct FieldDescriptor {
    std::string name, symbol;
    std::string domain="fluid", location="point";
    Parameters dimensions=Parameters::array({0,0,0,0,0,0,0});
    std::string storage="primary";
    std::string outputName;
    ValueType type=ValueType::Scalar;
    bool output=true, restart=false;
    Parameters initial, sets=Parameters::object(), boundaries=Parameters::object();
};
/// @brief 方程、约束、闭式、几何等对象的公共描述。
struct ObjectDescriptor {
    std::string name, category, type, expression;
    Parameters parameters=Parameters::object();
    /// @brief registry 的紧凑选择值，例如 multiPhase: eulerianEulerian。
    std::string selection;
};
/// @brief 文件、GUI、API 共用的中立输入。
struct Description {
    int formatVersion=1;
    std::string caseDir, name;
    std::vector<FieldDescriptor> fields;
    std::vector<ObjectDescriptor> objects;
    Parameters runtime=Parameters::object(), numerics=Parameters::object();
    Parameters solver=Parameters::object(), mesh=Parameters::object();
    /// @brief 按相对文件路径保存的文档注释，供格式 writer 原样带回。
    std::map<std::string,std::vector<std::string>> comments;
};
/// @brief 参数模式条目。模型注册时提供，IO 不枚举模型参数名。
struct ParameterRule {
    std::string name, type;
    bool required=false;
};
struct Schema {
    std::vector<ParameterRule> parameters;
    bool allowAdditional=false;
    void validate(const Parameters& value, const std::string& context) const;
};
/// @brief 通用、显式注册的工厂表；重复或未知类型直接报告。
template<class Context> class FactoryRegistry {
public:
    using Factory=std::function<void(const ObjectDescriptor&,Context&)>;
    void add(const std::string& type, Schema schema, Factory factory) {
        if (!entries_.emplace(type,Entry{std::move(schema),std::move(factory)}).second)
            throw std::runtime_error("Duplicate registered type: "+type);
    }
    bool contains(const std::string& type) const { return entries_.count(type)!=0; }
    const Schema& schema(const std::string& type) const { return entries_.at(type).schema; }
    void create(const ObjectDescriptor& object, Context& context) const {
        auto it=entries_.find(object.type);
        if(it==entries_.end()) throw std::runtime_error(object.category+" '"+object.name+"': unregistered type "+object.type);
        it->second.schema.validate(object.parameters,object.category+"/"+object.name);
        it->second.factory(object,context);
    }
private:
    struct Entry { Schema schema; Factory factory; };
    std::map<std::string,Entry> entries_;
};
/// @brief 字段运行时视图，可绑定任意现有标量/向量/张量存储。
struct FieldView {
    FieldDescriptor descriptor;
    std::size_t size=0;
    std::function<double(std::size_t,int)> read;
    std::function<void(std::size_t,int,double)> write;
};
class FieldRegistry {
public:
    /// @brief 注册非拥有字段；字段存储生命周期由模型保证。
    void add(FieldView field);
    FieldView& lookup(const std::string& name);
    const std::map<std::string,FieldView>& fields() const { return fields_; }
private:
    std::map<std::string,FieldView> fields_;
};
} // namespace SF::Model
