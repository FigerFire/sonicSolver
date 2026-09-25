#pragma once
#include "infrastructure/io/SF_caseConfig.h"
namespace SF::Application::ModelLoader {
/// @brief 把物理模型的场视图与通用输出管线绑定。
void configureOutput(ResultWriter& writer,const CaseConfig& config);
}
