#pragma once
#include <string>
namespace SF { struct MeshParameters;
namespace MeshIO {
/// @brief 读取结构网格生成参数，不依赖求解器配置。
bool readParameters(const std::string& path, MeshParameters& result);
}}
