/// @file SF_cli.h
/// @brief 命令行入口、参数解析与应用启动职责。

#pragma once

namespace SF::CLI {

/// @brief 命令行前端调度入口。
///
/// 支持 `run`、`check`、`explain`、`doctor`、`models`、`why`、
/// `explainModel`、`recipes`、`init`、`cleanResult`、`postProcessing` 和 `help`；不带
/// case 参数时统一使用当前目录 `.`，另提供 `initialOutput` 和
/// `steps` 两个运行控制命令。
/// @param argc 命令行参数数量。
/// @param argv 命令行参数。
/// @return 求解任务退出码。
int run(int argc, char* argv[]);

} // namespace SF::CLI
