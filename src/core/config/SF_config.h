#pragma once
/// @file SF_config.h
/// @brief config 目录总入口：保留历史兼容 include 面。
///
/// 配置类型按模块 ownership 拆分。字符串解析属于 application/model：core
/// 不能知道 application 存在，因此这里只重导出 core 的 typed 类型；调用方
/// 需要 parser 时显式 include "app/application/model/SF_configParser.h"。

#include "core/config/SF_configTypes.h"
#include "core/config/SF_numericsPolicy.h"
