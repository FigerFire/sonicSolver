#pragma once
#include "core/model/SF_model.h"
namespace SF::FieldIO {
/// @brief 通用字段描述序列化，字段名称无物理语义。
Model::FieldDescriptor decode(const std::string& name,const Model::Parameters& value);
/// @brief 解析 fields/registry 文件中的完整字段声明。
Model::FieldDescriptor decodeRegistry(const std::string& name,const Model::Parameters& value);
Model::Parameters encode(const Model::FieldDescriptor& field);
/// @brief 编码 fields/registry 文件中的字段声明，不包含内部场数值。
Model::Parameters encodeRegistry(const Model::FieldDescriptor& field);
/// @brief 按注册表读取/写入数值数组，组件数由基础类型确定。
void write(const Model::FieldRegistry& fields,const std::string& file,bool restartOnly=false);
void read(Model::FieldRegistry& fields,const std::string& file,bool restartOnly=false);
}
