#pragma once
#include <string>
namespace SF::Application {
/// @brief CLI/GUI 共用的原生工程创建；只支持真实实现的 recipe，拒绝非空目录。
void createNativeCase(const std::string& directory,const std::string& recipe="compressibleSingleFluid");
}
