#pragma once
#include "infrastructure/io/SF_caseConfig.h"
#include "core/model/SF_model.h"
namespace SF::Application::ModelLoader {
std::vector<EquationSystemInstanceConfig> equationInstances(
    const Model::Description& model);
/// @brief 模型/插件启动时注册自己的 schema 与强类型构造函数，不修改 CaseReader。
Model::FactoryRegistry<CaseConfig>& objectRegistry();
/// @brief 通用描述编译到现有内置求解器的强类型配置；文件 IO 不参与模型选择。
CaseConfig build(const Model::Description& model);
CaseConfig read(const std::string& path);
}
