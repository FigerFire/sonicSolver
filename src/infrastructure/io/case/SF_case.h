#pragma once
#include "core/model/SF_model.h"
namespace SF::CaseIO {
/// @brief 通过 object/type IORegistry 反序列化 SonicFile case，不创建任何物理模型。
Model::Description read(const std::string& path);
/// @brief 写出独立的规范 case 元数据；网格与大数组由各自 IO 管理。
void write(const Model::Description& model, const std::string& directory);
}
