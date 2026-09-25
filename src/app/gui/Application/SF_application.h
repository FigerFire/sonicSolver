/// @file SF_application.h
/// @brief GUI 应用生命周期与主窗口启动职责。

#pragma once

namespace SF::GUI::Application {

/// @brief 启动 sonicGui 应用。
/// @param argc 命令行参数数量。
/// @param argv 命令行参数。
/// @return GUI 事件循环退出码。
int run(int argc, char* argv[]);

} // namespace SF::GUI::Application
