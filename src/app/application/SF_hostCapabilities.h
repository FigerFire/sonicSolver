#pragma once

/// @file SF_hostCapabilities.h
/// @brief composition root 对 host backend 能力的唯一查询入口。

#include "core/interfaces/SF_buildCapabilities.h"

namespace SF::Application {

/// @brief 从真实链接/注册的 backend 推出本 binary 的能力。
///
/// 调用点只有 composition root：SystemBuilder 收到的是这份值的副本，
/// 不允许自己假定任何 backend 存在。
System::BuildCapabilities detectBuildCapabilities();

} // namespace SF::Application
