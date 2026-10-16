#pragma once
#include <string>
#include "core/model/SF_model.h"
namespace SF { struct MeshParameters;
namespace MeshIO {
/// @brief 读取结构网格生成参数，不依赖求解器配置。
void decodeParameters(const Model::Parameters& input, MeshParameters& result);
bool readParameters(const std::string& path, MeshParameters& result);
}}
