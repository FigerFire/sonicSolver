#pragma once
#include "core/model/SF_model.h"
namespace SF::Serialization {
/// @brief 已解析文档及其 YAML 注释。注释不参与模型语义校验。
struct Document {
    Model::Parameters value;
    std::vector<std::string> comments;
};
/// @brief SonicFile 的统一文件身份和正文。标题既可写在 SonicFile 节内，也兼容旧的顶层写法。
struct SonicDocument {
    std::string object, type, source;
    Model::Parameters body;
    std::vector<std::string> comments;
};
/// @brief YAML/JSON 与通用节点互转，保留 mapping 顺序并拒绝重复键。
Model::Parameters read(const std::string& path);
void write(const Model::Parameters& value,const std::string& path);
/// @brief 读取并收集行注释；注释会在写回时保留在文档头部。
Document readDocument(const std::string& path);
/// @brief 写入模型并保留 readDocument 收集的注释。
void writeDocument(const Document& document,const std::string& path);
/// @brief 读取带 object/type 标题的 SonicFile。支持紧凑的 condition: value 边界语法。
SonicDocument readSonicFile(const std::string& path);
/// @brief 写入可读的 SonicFile；compactConditions 将单键条件写为 patch: condition: value。
void writeSonicFile(const SonicDocument& document,const std::string& path,
                    bool compactConditions=false);
}
